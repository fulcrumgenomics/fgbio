/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.DownsampleAndNormalizeBam.CoverageTracker
import com.fulcrumgenomics.bam.api._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sorter}
import htsjdk.samtools.util.IntervalList.IntervalMergerIterator
import htsjdk.samtools.util.{CoordMath, Interval, IntervalList, Murmur3, OverlapDetector}

/**
  * Companion object with helper classes.
  */
private object DownsampleAndNormalizeBam extends LazyLogging {
  /**
    * Tracks coverage over a set of regions.
    */
  class CoverageTracker(regions: Seq[Interval],
                        val targetCoverage: Int,
                        val minMapQ: Int
                       ) {

    // Overlap detector for region => coverage
    private val detector = new OverlapDetector[RegionCoverage](0, 0)
    regions.foreach { region =>
      detector.addLhs(new RegionCoverage(region.getContig, region.getStart, region.getEnd), region)
    }

    /**
      * Checks to see if a template has reads that contribute coverage to any bases that are currently
      * below the coverage threshold.  If so, coverage from all acceptable reads is accumulated and the
      * method returns true.  Otherwise no coverage is accumulated and the method returns false.
      *
      * @param t a Template to accumulate coverage from if it covers any bases that are under the coverage threshold
      * @return true if the template added coverage, false otherwise
      */
    def add(t: Template): Boolean = {
      val recs = t.allReads
        .filterNot(_.secondary)
        .filterNot(_.unmapped)
        .filterNot(_.duplicate)
        .filter(_.mapq >= minMapQ)
        .toSeq

      val addsCoverage = recs.exists { rec =>
        detector.getOverlaps(rec.asSam).iterator().exists { cov =>
          Range.inclusive(math.max(rec.start, cov.start), math.min(rec.end, cov.end)).exists { pos =>
            cov(pos) <= targetCoverage
          }
        }
      }

      if (addsCoverage) {
        val intervals = {
          val iter = recs
            .sortBy(r => (r.refName, r.start, r.end))
            .iterator
            .map(r => new Interval(r.refName, r.start, r.end))
          new IntervalMergerIterator(iter, true, false, false)
        }

        for (interval <- intervals; cov <- detector.getOverlaps(interval).iterator()) {
          cov.add(interval.getStart, interval.getEnd)
        }
      }

      addsCoverage
    }
  }

  /**
    * Coverage counter for a region - all coordinates in the public API are 1-based closed ended.
    */
  class RegionCoverage(val chrom: String, val start: Int, val end: Int) {
    private val offset = start
    private val counts = new Array[Int](CoordMath.getLength(start, end))

    /** Returns the coverage at the given genomic position, or -1 if the position is not between start:end. */
    def apply(i: Int): Int = {
      if (i < start || i > end) -1
      else this.counts(i - offset)
    }

    /** Adds 1 to the coverage counts for all bases between from:to inclusive. */
    def add(from: Int, to: Int): Unit = {
      forloop (from=math.max(0, from-offset), until=math.min(counts.length, to-offset+1)) { i =>
        this.counts(i) += 1
      }
    }
  }

  /** Generates an iterator over the templates in the SamSource that is in a random order. */
  // TODO: move into Bams?
  def templateRandomIterator(in: SamSource,
                             randomSeed: Int = 42,
                             maxInMemory: Int = Bams.MaxInMemory,
                             tmpDir: DirPath = Io.tmpDir): Iterator[Template] = {
    val hasher = new Murmur3(randomSeed)
    val sorter = {
      val f = (r: SamRecord) => SamOrder.RandomQueryKey(hasher.hashUnencodedChars(r.name), r.name, r.asSam.getFlags)
      new Sorter(maxInMemory, new SamRecordCodec(in.header), f)
    }

    val sortProgress = ProgressLogger(logger, verb = "sorted", unit = 5e6.toInt)
    in.foreach { rec =>
      sorter += rec
      sortProgress.record()
    }

    val header = in.header.clone()
    SamOrder.RandomQuery.applyTo(header)
    Bams.templateIterator(sorter.iterator, header, maxInMemory, tmpDir)
  }
}

@clp(
  group=ClpGroups.SamOrBam,
  description =
    """
      |Downsamples a BAM in a biased way to a uniform coverage across regions.
      |
      |Attempts to downsample a BAM such that every base in the genome (or in the target `regions` if provided)
      |is covered by at least `coverage` reads.  When computing coverage:
      |  - Reads marked as secondary, duplicate or unmapped are not used
      |  - A base can receive coverage from only one read with the same queryname (i.e. mate overlaps are not counted)
      |  - Coverage is counted if a read _spans_ a base, even if that base is deleted in the read
      |
      |Reads are first sorted into a random order (by hashing read names).  Reads are then consumed one template
      |at a time, and if any read adds coverage to base that is under the target coverage, _all_ reads (including
      |secondary, unmapped, etc.) for that template are emitted into the output.
      |
      |Given the procedure used for downsampling, it is likely the output BAM will have coverage up to 2X the requested
      |coverage at regions in the input BAM that are i) well covered and ii) are close to regions that are poorly
      |covered.
    """
)
class DownsampleAndNormalizeBam
( @arg(flag='i', doc="Input SAM or BAM file.") val input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam,
  @arg(flag='c', doc="Desired minimum coverage.") val coverage: Int,
  @arg(flag='m', doc="Minimum mapping quality to count a read as covering.") val minMapQ: Int = 0,
  @arg(flag='s', doc="Random seed to use when randomizing order of reads/templates.") val seed: Int = 42,
  @arg(flag='l', doc="Optional set of regions for coverage targeting.") val regions: Option[PathToIntervals] = None,
  @arg(flag='M', doc="Maximum records to be held in memory while sorting.") val maxInMemory: Int = Bams.MaxInMemory
) extends FgBioTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)
  regions.foreach(Io.assertReadable)

  override def execute(): Unit = {
    val in  = SamSource(input)
    val out = SamWriter(output, in.header, sort=SamOrder(in.header))

    val targets: IndexedSeq[Interval] = regions match {
      case Some(rs) => IntervalList.fromPath(rs).uniqued().getIntervals.toIndexedSeq
      case None     => in.dict.map(s => new Interval(s.name, 1, s.length)).toIndexedSeq
    }

    val tracker = new CoverageTracker(targets, this.coverage, this.minMapQ)
    var (read, written) = (0L, 0L)

    DownsampleAndNormalizeBam.templateRandomIterator(in, this.seed, maxInMemory).foreach { template =>
      read += 1

      if (tracker.add(template)) {
        out ++= template.allReads
        written += 1
      }

      if (read % 1000000 == 0) {
        logger.info(f"Read ${read}%,d templates and wrote out ${written}%,d")
      }
    }

    in.safelyClose()
    out.close()
    logger.info(f"Read ${read}%,d templates and wrote out ${written}%,d")
  }
}
