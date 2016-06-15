/**
 * Copyright (c) 2015, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

package com.fulcrumgenomics.util;

import htsjdk.samtools.util.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import static java.lang.Math.pow;

/**
 * Program for picking sets of indices of arbitrary length that meet certain constraints
 * and attempt to maximize the edit distance between all members of the set picked.
 *
 * @author Tim Fennell
 */
class PickIlluminaIndicesCommand {
    /** The temperature at which the adapter ligation happens. */
    public static final int LIGATION_TEMPERATURE = 30;

    public int INDEX_LENGTH;
    public int N_INDICES;
    public int EDIT_DISTANCE;
    public boolean ALLOW_REVERSES;
    public boolean ALLOW_REVERSE_COMPLEMENTS;
    public boolean ALLOW_PALINDROMES;
    public int MAX_HOMOPOLYMER;
    public double MIN_GC;
    public double MAX_GC;
    public File OUTPUT;
    public int NUM_THREADS;
    public File VIENNA_RNA_DIR;
    public double MIN_DELTAG;
    public List<String> INDEX_ADAPTER;
    public List<String> AVOID_SEQUENCE;

    final static Log log = Log.getInstance(PickIlluminaIndicesCommand.class);

    private static final int[] SCORE_BY_DISTANCE = {0, 500000, 5000, 1};
    private int defaultEdgeSetSize = 0;

    /**
     *  Little class to encapsulate a index sequence and it's relationships to other indices.
     */
    class Index {
        final byte[] sequence;
        private final Collection<Index> related;
        private int _score;
        private byte _minEditDistance;
        private List<DnaFoldPrediction> folds = new ArrayList<>();

        Index(final byte[] sequence) {
            this.sequence = sequence;
            this.related = new HashSet<>(defaultEdgeSetSize);
            reset();
        }

        private void reset() { _score = Integer.MIN_VALUE; _minEditDistance = -1; }

        synchronized void addRelated(final Index b) {
            this.related.add(b);
            reset();
        }

        void removeRelated(final Index b) {
            this.related.remove(b);
            reset();
        }

        /** Returns the barcode that is this one reversed. */
        Index reverse() {
            final byte[] tmp = Arrays.copyOf(this.sequence, this.sequence.length);
            SequenceUtil.reverseQualities(tmp); // tmp isn't a quality score array, but all this method does is reverse in place!
            return new Index(tmp);
        }

        /** Generates a score where large scores are bad. */
        final int score() {
            if (_score == Integer.MIN_VALUE ) recalculate();
            return _score;
        }

        /** Returns the minimum edit distance between this barcode and any barcode in it's related set. */
        final int minEditDistance() {
            if (_minEditDistance == -1) recalculate();
            return _minEditDistance;
        }

        /** Gives the edit distance between this barcode and another of the same length. */
        byte calculateEditDistance(final Index that) {
            byte tmp = 0;
            for (int i=0; i<this.sequence.length; ++i) {
                if (this.sequence[i] != that.sequence[i]) ++tmp;
            }
            return tmp;
        }

        /** Recalculates the score and minEditDistance after changes are made. */
        final void recalculate() {
            _score = Integer.MIN_VALUE + 1;
            _minEditDistance = (byte) INDEX_LENGTH;
            for (final Index that : this.related) {
                final byte distance = calculateEditDistance(that);

                if (distance < SCORE_BY_DISTANCE.length) _score += SCORE_BY_DISTANCE[distance];
                if (distance < _minEditDistance) _minEditDistance = distance;
            }
        }

        final int getCountAtEditDistance(final int n) {
            int count = 0;
            for (final Index that : this.related) {
                final int distance = calculateEditDistance(that);
                if (distance == n) ++count;
            }

            return count;
        }

        /** Returns the barcode that is this one reverse complemented. */
        Index reverseComplement() {
            final byte[] tmp = Arrays.copyOf(sequence, sequence.length);
            SequenceUtil.reverseComplement(tmp);
            return new Index(tmp);
        }

        @Override public int hashCode() { return Arrays.hashCode(sequence); }

        @Override public boolean equals(final Object o) {
            return o != null && o instanceof Index && Arrays.equals(this.sequence, ((Index) o).sequence);
        }

        @Override public String toString() {
            return "{sequence: " + StringUtil.bytesToString(this.sequence) + ", score: " + score() + ", minEditDistance:" + minEditDistance() + "}";
        }
    }
    ///////////////////////////////////////////////////////////////////////////
    // End of Barcode class
    ///////////////////////////////////////////////////////////////////////////

    /**
     * A concurrent, multi-threaded, sorted set that does concurrent removal and updating of the ordering when a
     * barcode is removed.
     */
    class IndexSet extends ConcurrentSkipListSet<Index> {
        final BlockingQueue<Runnable> queue = new ArrayBlockingQueue<>(10000);
        final ThreadPoolExecutor exec = new ThreadPoolExecutor(NUM_THREADS, NUM_THREADS, 1, TimeUnit.MINUTES, queue);

        /** Basic constructor. */
        IndexSet(final Comparator<? super Index> comparator) {
            super(comparator);
        }

        /**
         * Overridden to remove a Barcode and the "edit" all the related indices to remove their relationship
         * to this barcode and re-insert them into the Set to resort them.
         */
        @Override public boolean remove(final Object o) {
            final Index index = (Index) o;

            /** Initialize a countdown int that each thread will update so we know when everything is processed. */
            final AtomicInteger countdown = new AtomicInteger(index.related.size());
            final Object latch = new Object();

            if (super.remove(index)) {
                /** For each barcode queue up a job to recalculate scores and remove/add from the ranked set. */
                for (final Index other : index.related) {
                    exec.execute(() -> {
                        IndexSet.super.remove(other);
                        other.removeRelated(index);
                        other.recalculate();
                        IndexSet.this.add(other);

                        if (countdown.decrementAndGet() == 0) {
                            synchronized (latch) { latch.notify(); }
                        }
                    });
                }

                /** Wait until all the related indices get updated. */
                while (countdown.get() > 0) {
                    try { synchronized(latch) { latch.wait(1000); } }
                    catch (InterruptedException ie) { /* Do Nothing */ }
                }

                return true;
            }
            else {
                return false;
            }
        }
    }

    protected int execute() {
        // Generate all kmers
        final List<byte[]> kmers = generateAllKmers(INDEX_LENGTH);
        List<Index> indexes = new LinkedList<>();
        log.info("Generated " + kmers.size() + " kmers.");

        this.defaultEdgeSetSize = guessAtNumberOfEdgesPerBarcode(INDEX_LENGTH, EDIT_DISTANCE);
        log.info("Guess at number of edges: " + this.defaultEdgeSetSize);

        // Remove any that violate the max homopolyer length and convert them into indices
        for (final byte[] kmer : kmers) {
            if (lengthOfLongestHomopolymer(kmer) <= MAX_HOMOPOLYMER) {
                indexes.add(new Index(kmer));
            }
        }
        log.info("Retained " + indexes.size() + " indices after applying max homopolymer restriction.");

        if (!ALLOW_PALINDROMES) {
            filterForPalindromes(indexes);
            log.info("Retained " + indexes.size() + " indices after filtering for palindromes.");
        }

        filterForAvoidSequences(indexes);
        log.info("Retained " + indexes.size() + " indices after filtering for sequences to avoid.");

        filterForGc(indexes);
        log.info("Retained " + indexes.size() + " indices after filtering for GC content.");

        // Filter out reverses
        if (!ALLOW_REVERSES) {
            filterOutReverses(indexes);
            log.info("Retained " + indexes.size() + " indices after removing indices that are reversals of other indices.");
        }

        // Filter out reverse complements
        if (!ALLOW_REVERSE_COMPLEMENTS) {
            filterOutReverseComplements(indexes);
            log.info("Retained " + indexes.size() + " indices after removing indices that are reverse-complements of other indices.");
        }

        // Filter out indices that have very bad secondary structure
        {
            DnaFoldPredictor predictor = new DnaFoldPredictor(VIENNA_RNA_DIR, LIGATION_TEMPERATURE);
            final Iterator<Index> iterator = indexes.iterator();
            while (iterator.hasNext()) {
                final Index index = iterator.next();
                final String bases = StringUtil.bytesToString(index.sequence);

                // Predict folding for each adapter that the index will be embedded in
                for (final String adapter : this.INDEX_ADAPTER) {
                    final String seq = adapter.replaceFirst("N+", bases);
                    final DnaFoldPrediction prediction = predictor.predict(seq);
                    index.folds.add(prediction);
                }

                final double deltaG = index.folds.stream().mapToDouble(DnaFoldPrediction::deltaG).min().orElse(MIN_DELTAG);
                if (deltaG < MIN_DELTAG) iterator.remove();
            }
            log.info("Retained " + indexes.size() + " indices after removing indicies with too low deltaG.");
        }

        if (indexes.isEmpty()) {
            log.error("No indices left after previous filters.");
            System.exit(1);
        }

        // Make a big graph that connects all indices with edit distances between 1 and the
        // maximum EDIT_DISTANCE + 1
        log.info("Building the graph of all indices.");
        buildGraph(indexes);

        // First build a tree set where the "first" item is the one with the highest number of other indices
        log.info("Ranking indices.");
        final IndexSet rankedIndexes = rankBarcodes(indexes);
        indexes = null;

        // Finally go through and throw out indices in order until we achieve the desired edit distance
        int lastMinEdit = 0;
        int count = rankedIndexes.size();
        while (count > N_INDICES || lastMinEdit < EDIT_DISTANCE) {
            final Index first = rankedIndexes.first();

            if (first.minEditDistance() > lastMinEdit) {
                lastMinEdit = first.minEditDistance();
                log.info("There are " + count + " indices with a minimum edit distance of " + lastMinEdit);
            }

            rankedIndexes.remove(first);
            if (--count % 100 == 0) log.info("Down to " + count + " indices.");
        }

        log.info("Ended with " + count + " indices with a min edit distance of " + rankedIndexes.first().minEditDistance());
        try {
            final int adapters = this.INDEX_ADAPTER.size();
            final BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
            final NumberFormat fmt = NumberFormat.getIntegerInstance();
            final NumberFormat pct = NumberFormat.getPercentInstance();
            final NumberFormat dec = new DecimalFormat("0.00");
            out.write("Barcode\tMinEditDistance\tNumBarcodesAtMinEditDistance\tGC%");
            for (int i=1; i<=adapters; ++i) {
                for (final String header : new String[] {"seq", "fold", "deltag"}) {
                    out.write("\t");
                    out.write(header + "_" + i);
                }
            }
            out.newLine();

            for (final Index bc : rankedIndexes) {
                final int minEditDistance = bc.minEditDistance();
                final int numberAtEditDistance = bc.getCountAtEditDistance(minEditDistance);

                out.write(StringUtil.bytesToString(bc.sequence));
                out.write('\t');
                out.write(fmt.format(minEditDistance));
                out.write('\t');
                out.write(fmt.format(numberAtEditDistance));
                out.write('\t');
                out.write(pct.format(SequenceUtil.calculateGc(bc.sequence)));
                for (int i=0; i<adapters; ++i) {
                    out.write('\t');
                    out.write(bc.folds.get(i).sequence());
                    out.write('\t');
                    out.write(bc.folds.get(i).structure());
                    out.write('\t');
                    out.write(dec.format(bc.folds.get(i).deltaG()));
                }
                out.newLine();
            }

            out.flush();
            out.close();
        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }

        return 0;
    }

    /**
     * Takes a List of indices and removes any that are outside the desired GC range.
     */
    private void filterForGc(final List<Index> indexes) {
        final Iterator<Index> iterator = indexes.iterator();
        while (iterator.hasNext()) {
            final Index b = iterator.next();
            final double gc = SequenceUtil.calculateGc(b.sequence);
            if (gc < MIN_GC || gc > MAX_GC) iterator.remove();
        }
    }

    /**
     * Takes the list of indices and identifies ones that have reverse complements in the set and then removes
     * from the set the barcode that comes later in the list.
     */
    private void filterOutReverseComplements(final List<Index> indexes) {
        final LinkedHashSet<Index> tmp = new LinkedHashSet<>();
        for (final Index bc : indexes) {
            if (!tmp.contains(bc.reverseComplement())) tmp.add(bc);
        }
        indexes.clear();
        indexes.addAll(tmp);
    }

    /**
     * Takes the list of indices and identifies ones that have reverses (not RCs) in the set and then removes
     * from the set the barcode that comes later in the list.
     */
    private void filterOutReverses(final List<Index> indexes) {
        final LinkedHashSet<Index> tmp = new LinkedHashSet<>();
        for (final Index bc : indexes) {
            if (!tmp.contains(bc.reverse())) tmp.add(bc);
        }
        indexes.clear();
        indexes.addAll(tmp);
    }

    /** Filters out any indices that are kmers contained in the list of sequences to avoid. */
    private void filterForAvoidSequences(final List<Index> indexes) {
        final Set<String> avoidKmers = new HashSet<>();
        for (final String seq : AVOID_SEQUENCE) {
            for (int i=0; i<seq.length() - INDEX_LENGTH + 1; ++i) {
                final String sub = seq.substring(i, i+ INDEX_LENGTH);
                avoidKmers.add(sub);
                avoidKmers.add(SequenceUtil.reverseComplement(sub));
            }
        }

        final Iterator<Index> iterator = indexes.iterator();
        while (iterator.hasNext()) {
            final String barcodeString = StringUtil.bytesToString(iterator.next().sequence);
            if (avoidKmers.contains(barcodeString)) {
                iterator.remove();
            }
        }
    }

    /** Removes all sequences where the sequence == revcomp(sequence) */
    private void filterForPalindromes(final List<Index> indexes) {
        final Iterator<Index> iterator = indexes.iterator();
        while (iterator.hasNext()) {
            final Index bc = iterator.next();
            final byte[] rc = Arrays.copyOf(bc.sequence, bc.sequence.length);
            SequenceUtil.reverseComplement(rc);
            if (Arrays.equals(bc.sequence, rc)) {
                iterator.remove();
            }
        }
    }

    /**
     * Generates a SortetSet of indices where they are ranked by the "most connected" to
     * the "least connected".  Must be called after constructing the graph.
     */
    private IndexSet rankBarcodes(final Collection<Index> input) {
        final Comparator<Index> comparator = (lhs, rhs) -> {
            // Remember: BIGGER scores == WORSE scores
            int retval = rhs.score() - lhs.score();
            if (retval == 0) {
                for (int i=0; i<lhs.sequence.length && retval == 0; ++i) retval = lhs.sequence[i] - rhs.sequence[i];
            }

            return retval;
        };

        final IndexSet set = new IndexSet(comparator);
        set.addAll(input);
        return set;
    }

    /**
     * Builds a graph between all the indices in the list provided.  The indices are connected in the graph if they
     * have EDIT_DISTANCE or fewer edits between the pair.  The graph is constructed using a thread pool so that
     * it can run in a reasonable time.
     */
    private void buildGraph(final List<Index> indexes) {
        final ExecutorService exec = Executors.newFixedThreadPool(NUM_THREADS);

        final Index[] bcs = indexes.toArray(new Index[indexes.size()]);
        final int len = bcs.length;
        final AtomicInteger counter = new AtomicInteger(0);

        for (int i=0; i<len; ++i) {
            final Index lhs = bcs[i];
            final int firstJ = i+1;

            exec.submit((Runnable) () -> {
                for (int j=firstJ; j<len; ++j) {
                    final Index rhs = bcs[j];
                    final int distance = lhs.calculateEditDistance(rhs);
                    if (distance <= EDIT_DISTANCE) {
                        lhs.addRelated(rhs);
                        rhs.addRelated(lhs);
                    }
                }

                final int count = counter.incrementAndGet();
                if (count % 1000 == 0) {
                    log.info("    Processed " + count + " indices.");
                }
            });
        }

        try {
            exec.shutdown();
            exec.awaitTermination(7, TimeUnit.DAYS);
        }
        catch (InterruptedException ie) {
            throw new RuntimeException(ie);
        }
    }

    /** Generates all possible kmers of length and returns them as byte[]s. */
    List<byte[]> generateAllKmers(final int length) {
        final List<byte[]> sofar = new LinkedList<>();
        final byte[] bases = {'A', 'C', 'G', 'T'};

        if (sofar.size() == 0) {
            sofar.add(new byte[length]);
        }

        while (true) {
            final byte[] bs = sofar.remove(0);
            final int indexOfNextBase = findIndexOfNextBase(bs);

            if (indexOfNextBase == -1) {
                sofar.add(bs);
                break;
            }
            else {
                for (final byte b : bases) {
                    final byte[] next = Arrays.copyOf(bs, bs.length);
                    next[indexOfNextBase] = b;
                    sofar.add(next);
                }
            }
        }

        return sofar;
    }

    /** Finds the first non-zero character in the array, or returns -1 if all are non-zero. */
    int findIndexOfNextBase(final byte[] bs) {
        for (int i=0; i<bs.length; ++i) {
            if (bs[i] == 0) return i;
        }

        return -1;
    }

    // Finds the longest homopolymer in a stretch of bases.  Not terribly efficient!
    int lengthOfLongestHomopolymer(final byte[] bytes) {
        int longest = 1;

        for (int i=0; i<bytes.length; ++i) {
            int length = 0;

            for (int j=i; j<bytes.length; ++j) {
                if (bytes[i] == bytes[j]) {
                    ++length;
                }
                else {
                    break;
                }
            }

            longest = Math.max(longest, length);
        }

        return longest;
    }

    /** Guesses at how many connections each barcode will have, in order to presize a lot of sets. */
    int guessAtNumberOfEdgesPerBarcode(final int length, final int edits) {
        int total = 0;
        for (int i=edits; i>0; --i) total += nPickK(length, i) * (int) pow(3, i);
        return (int) (total / 2d);
    }

    /** Implementation of the calculation for how many combinations are there of K items from N items. */
    final int nPickK(final int n, final int k) {
        return fac(n) / ( fac(k)  * fac(n-k));
    }

    /** Simple factorial() implementation. */
    final int fac(final int n) {
        int result = n;
        for (int i=n-1; i>1; --i) result *= i;
        return result;
    }
}

