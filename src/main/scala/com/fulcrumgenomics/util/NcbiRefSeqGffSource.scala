/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

import java.lang.Integer.{parseInt => int}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.util.GeneAnnotations.{Exon, Gene, GeneLocus, Transcript}
import htsjdk.samtools.{SAMSequenceDictionary, SAMSequenceRecord}

import scala.collection.mutable

/**
  * Companion object for [[NcbiRefSeqGffSource]] which provides factory methods for parsing RefSeq GFF
  * files. The resulting data is returned as instances of [[GeneAnnotations.Gene]] and related classes.
  *
  * The following link points to a file of the appropriate format that is parsed by this code:
  *   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Macaca_fascicularis/latest_assembly_versions/GCF_000364345.1_Macaca_fascicularis_5.0/GCF_000364345.1_Macaca_fascicularis_5.0_genomic.gff.gz
  */
object NcbiRefSeqGffSource {
  /** Represents a generic GFF feature with it's fixed set of fields (minus score).
    *
    * @param accession the accession of the sequence on which the feature is found
    * @param source the source of the sequence (e.g. RefSeq, Gnomon, etc.)
    * @param kind the kind of feature (e.g. gene, mRNA, exon, etc.)
    * @param start the start position of the feature (1-based)
    * @param end the end position of the feature (1-based inclusive)
    * @param negative whether the feature is on the negative strand
    * @param frame the coding frame of the feature. If the feature does not have a frame, is set to -1.
    * @param attrs the string of `;` separated `key=value` dynamic attributes on the record
    */
  private case class GffRecord(accession: String, source: String, kind: String, start: Int, end: Int, negative: Boolean, frame: Int, attrs: String) {
    // Lazily break the attributes out into a Map; lazy in case it's a feature type we don't need this for
    private lazy val attributes = attrs.split(';').map(_.split('=')).map { case Array(k,v) => k -> v }.toMap

    /** Gets a dynamic attribute by name. Throws an exception if the attribute is not present. */
    def apply(name: String): String = this.attributes(name)

    /** Returns `true` if the dynamic attribute exists on this record, `false` otherwise. */
    def has(name: String): Boolean = this.attributes.contains(name)

    /** Returns the ID of this record from it's dynamic attributes. The ID is used to link records together. */
    def id: String = this("ID")

    /** Returns the ID of the parent record for this record. This is not populated for records of kind `gene`
      * but is for all sub-features of genes. Genes parent transcripts and transcript-like entries, transcripts
      * parent exons and CDS records.
      */
    def parentId: Option[String] = this.attributes.get("Parent")

    /** Returns the name of this record. All record types we parse have names, but not every GFF record has a name. */
    def name: String = this("Name")
  }

  /** Parses a RefSeq GFF file from an iterator of lines of text.
    *
    * @param lines the lines to parse
    * @param includeXs whether to include experimental transcripts (i.e. XM_* XP_* and XR_*).
    * @param dict a sequence dictionary used to help resolve accessions used in the GFF
    * @return an instance of [[NcbiRefSeqGffSource]] that can be queried for contents
    */
  def apply(lines: Iterator[String], includeXs: Boolean, dict: SAMSequenceDictionary): NcbiRefSeqGffSource =
    new NcbiRefSeqGffSource(lines, includeXs, dict)

  /** Parses a RefSeq GFF file from a collection of lines of text.
    *
    * @param lines the lines to parse
    * @param includeXs whether to include experimental transcripts (i.e. XM_* XP_* and XR_*).
    * @param dict a sequence dictionary used to help resolve accessions used in the GFF
    * @return an instance of [[NcbiRefSeqGffSource]] that can be queried for contents
    */
  def apply(lines: Iterable[String], includeXs: Boolean, dict: SAMSequenceDictionary): NcbiRefSeqGffSource =
    NcbiRefSeqGffSource(lines.iterator, includeXs, dict)

  /** Parses a RefSeq GFF file from an (optionally gzipped) text file.
    *
    * @param path the path to the file to parse
    * @param includeXs whether to include experimental transcripts (i.e. XM_* XP_* and XR_*).
    * @param dict a sequence dictionary used to help resolve accessions used in the GFF
    * @return an instance of [[NcbiRefSeqGffSource]] that can be queried for contents
    */
  def apply(path: FilePath, includeXs: Boolean, dict: SAMSequenceDictionary): NcbiRefSeqGffSource =
    NcbiRefSeqGffSource(Io.readLines(path), includeXs, dict)
}


/** Queryable object that contains gene/transcript/exon information parsed from a RefSeq GFF file.
  *
  * One notable complication in parsing RefSeq GFF files is that they refer to chromosomes/contigs by NCBI
  * accessions rather than by more common names (e.g. chr1).  The accessions are translated such that the returned
  * records use common contig names.  The `dict` parameter is the preferred way to do this.  Ideally a dictionary
  * should be provided with accessions as sequence aliases, e.g.:
  *
  *   @SQ     SN:chr1 LN:248956422    M5:2648ae1bacce4ec4b6cf337dcae37816     AS:hg38 AN:NC_000001.11
  *
  * This will allow the chromosome to be looked up by accession and resolved to it's common name without any
  * other assumptions.  If an accession in the GFF is not found in the sequence dictionary it will be handled as
  * follows:
  *
  *   - NC_* accessions with name MT in the GFF are translated to chrM
  *   - NC_* accessions other than MT are automatically mapped to "chr{chromsome-number}"
  *   - Other accessions/contigs are skipped
  *
  * @param lines the GFF lines to parse
  * @param includeXs whether to include experimental transcripts (i.e. XM_* XP_* and XR_*).
  * @param dict a sequence dictionary used to help resolve accessions in the accessions in the GFF
  */
class NcbiRefSeqGffSource private(lines: Iterator[String],
                                  val includeXs: Boolean,
                                  val dict: SAMSequenceDictionary) extends Iterable[Gene] with LazyLogging {
  import NcbiRefSeqGffSource._

  // Construct a lookup from the sequence dictionary that includes all the aliases as well as primary names
  private val sequenceLookup: Map[String, SAMSequenceRecord] = dict.getSequences.iterator().flatMap { s =>
    s.getAttribute("AN") match {
      case null    => Seq(s.getSequenceName -> s)
      case aliases => (s.getSequenceName +: aliases.split(',').toSeq).map(name => name -> s)
    }
  }.toMap

  private val genes  = mutable.LinkedHashMap[String,Gene]()

  /** NOTE about the structure of the RefSeq GFF files.
    *
    * Although GFF files are inherently flat with each line representing a record, it is best to
    * think of RefSeq GFF files as hierarchical files with different kinds of records "nested"
    * within each other.  The hierarchy looks something like:
    *   * Chromosome region line
    *       * Gene lines
    *           * Transcript (and transcript-like entity) lines
    *               * Exon lines
    *               * CDS lines
    *
    * Thus the parsing takes place top-down by a) scanning for the next line of the expected type
    * and then b) attempting to parse that line and it's children.
    */
  private val _bufferedLines = lines.bufferBetter

  /** The genome build declared in the header of the GFF. */
  val genomeBuild: Option[String] = _bufferedLines.find(_.startsWith("#!genome-build")).map(_.split(' ')(1))

  {
    val iter = _bufferedLines.filter(l => !l.startsWith("#")).map { line =>
      val fields = line.split('\t')
      GffRecord(
        accession = fields(0),
        source    = fields(1),
        kind      = fields(2),
        start     = int(fields(3)),
        end       = int(fields(4)),
        negative  = fields(6) == "-",
        frame     = if (fields(7) == ".") -1 else int(fields(7)),
        attrs     = fields(8))
    }.bufferBetter

    while (iter.nonEmpty) {
      // Scan for the next chromosome record, discard any other records along the way
      iter.find(rec => rec.kind == "region") match {
        case None =>
          () // Indicates we've hit the end of the file
        case Some(rec) =>
          val chromName = this.sequenceLookup.get(rec.accession) match {
            case Some(sequence)                          => Some(sequence.getSequenceName)
            case None if rec.name == "MT"                => Some("chrM")
            case None if rec.accession.startsWith("NC_") => Some("chr" + rec.name)
            case None                                    => None
          }

          chromName match {
            case None =>
              logger.debug(s"Skipping records for reference sequence accession ${rec.accession}")
              iter.dropWhile(_.accession == rec.accession)
            case Some(chrom) =>
              // Build a sub-iterator that stops at the end of this chromosome
              val chromIter = iter.takeWhile(_.accession == rec.accession)

              // Run through the chromosome's records, parsing out genes
              chromIter.foreach { rec =>
                if (rec.kind == "gene") {
                  try {
                    parseGene(chrom, rec, iter).foreach { gene =>
                      this.genes.get(gene.name) match {
                        case None =>
                          this.genes.put(gene.name, gene)
                        case Some(existing) =>
                          logger.debug(s"Adding second locus for $gene.")
                          this.genes.put(gene.name, existing.copy(loci = existing.loci ++ gene.loci))
                      }
                    }
                  }
                  catch {
                    case ex: Exception => logger.error(s"Error parsing gene: $rec"); throw ex
                  }
                }
              }
          }
      }
    }
  }
  logger.debug(s"Loaded ${genes.size} genes.")

  /**
    * Parses a single gene and all of it's sub-features from the iterator. After invoking this method
    * the `iter` passed in will be positioned such that the next record is the record _after_ all
    * children of the gene record (and their children and so on).
    *
    * @return a [[Some(Gene)]] if one of more transcripts are parsed, or [[None]] if no transcripts are
    *         parsed. A [[None]] is primarily returned if the gene has only experimentally determined
    *         transcript accessions (X?_*), and [[includeXs]] is false.
    * */
  private def parseGene(chrom: String, geneRec: GffRecord, iter: BetterBufferedIterator[GffRecord]): Option[Gene] = {
    val txsBuilder = Seq.newBuilder[Transcript]

    // Loop while the next record is a child of this gene and is "transcript-like", which includes any
    // direct children that define the attribute "transcript_id", not just records with kind=transcript.
    while (iter.hasNext && iter.head.parentId.contains(geneRec.id) && iter.head.has("transcript_id")) {
      val txRec = iter.next()
      val tx = parseTranscript(chrom, txRec, iter)
      if (this.includeXs || !tx.name.startsWith("X")) txsBuilder += tx
    }

    val txs = txsBuilder.result()
    if (txs.isEmpty) None else {
      val locus = GeneLocus(txs)
      Some(Gene(geneRec.name, Seq(locus)))
    }
  }

  /** Parses out a transcript and all of it's child records. After invoking this method the `iter` passed
    * in will be positioned such that the next records is immediately after the last child of this
    * transcript.
    */
  private def parseTranscript(chrom: String, txRec: GffRecord, iter: BetterBufferedIterator[GffRecord]): Transcript = {
    // Consume _all_ the children, including ones we may want to ignore
    val features = iter.takeWhile(_.parentId.contains(txRec.id)).toIndexedSeq

    // Build up the exons of the transcript and sort them into transcription order
    val exons = features.filter(_.kind == "exon").map(rec => Exon(rec.start, rec.end)).sortBy(e => if (txRec.negative) -e.end else e.start)

    // Find the coding start and end from the CDS features
    val (cdsStart, cdsEnd) = features.filter(_.kind == "CDS") match {
      case Seq() => (None, None)
      case xs    => (Some(xs.minBy(_.start).start), Some(xs.maxBy(_.end).end))
    }

    Transcript(txRec("transcript_id"), chrom, txRec.start, txRec.end, cdsStart, cdsEnd, txRec.negative, exons)
  }

  /**
    * Looks up a gene by name. Only recognizes the gene names from the GFF without any aliasing.
    * @param name the name of the gene
    * @return either Some(gene) if the gene exists or None otherwise.
    */
  def get(name: String): Option[Gene] = this.genes.get(name)

  /** Provides an iterator over the genes within the source.  Iteration is ordered by the first locus of each gene
    * but generally the results should be sorted as desired if sort order is important.
    */
  def iterator: Iterator[Gene] = this.genes.values.iterator
}
