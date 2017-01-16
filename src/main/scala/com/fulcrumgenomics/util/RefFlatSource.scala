/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics
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
 *
 */

package com.fulcrumgenomics.util

import java.io.{Closeable, File, InputStream}

import com.fulcrumgenomics.util.GeneAnnotations.{Exon, Gene, Transcript}
import dagr.commons.CommonsDef.FilePath
import dagr.commons.util.LazyLogging
import htsjdk.samtools.SAMSequenceDictionary

import scala.io.Source

object RefFlatSource {
  /** Creates a new refflat source from a sequence of lines. */
  def apply(lines: Seq[String], dict: Option[SAMSequenceDictionary]): RefFlatSource = new RefFlatSource(lines.iterator, dict)

  /** Creates a new fasrefflattq source from an iterator of lines. */
  def apply(lines: Iterator[String], dict: Option[SAMSequenceDictionary]): RefFlatSource = new RefFlatSource(lines, dict)

  /** Creates a new refflat source from an input stream. */
  def apply(stream: InputStream, dict: Option[SAMSequenceDictionary]): RefFlatSource = new RefFlatSource(Source.fromInputStream(stream).getLines(), dict, Some(stream))

  /** Creates a new refflat source from a source. */
  def apply(source: Source, dict: Option[SAMSequenceDictionary]): RefFlatSource = new RefFlatSource(source.getLines(), dict, Some(source))

  /** Creates a new refflat source from a File. */
  def apply(file: File, dict: Option[SAMSequenceDictionary]): RefFlatSource = apply(path=file.toPath, dict=dict)

  /** Creates a new refflat source from a Path. */
  def apply(path: FilePath, dict: Option[SAMSequenceDictionary] = None): RefFlatSource = apply(Io.toInputStream(path), dict)

  private val Header: Seq[String] = Seq(
    "geneName",
    "name",
    "chrom",
    "strand",
    "txStart",
    "txEnd",
    "cdsStart",
    "cdsEnd",
    "exonCount",
    "exonStarts",
    "exonEnds"
  )

  private val PicardHeader: Seq[String] = Seq(
    "GENE_NAME",
    "TRANSCRIPT_NAME",
    "CHROMOSOME",
    "STRAND",
    "TX_START",
    "TX_END",
    "CDS_START",
    "CDS_END",
    "EXON_COUNT",
    "EXON_STARTS",
    "EXON_ENDS"
  )

  private val PicardHeaderMap: Map[String, String] = PicardHeader.zip(Header).toMap
}

/** Reads gene annotation information from a RefFlat file.
  *
  * Skips genes on unrecognized chromosomes if a sequence dictionary is provided.
  *
  * The format is described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat
  *
  * A Picard-style header is also supported (GENE_NAME, TRANSCRIPT_NAME, ...).
  * */
class RefFlatSource private(lines: Iterator[String],
                            dict: Option[SAMSequenceDictionary] = None,
                            private[this] val source: Option[{ def close(): Unit }] = None
                           ) extends Iterable[Gene] with Closeable with LazyLogging {
  private val genes = {

    // Ensure there is a header for DelimitedDataParser, and if a Picard-style header is used, convert it
    val _lines: Iterator[String] = {
      import com.fulcrumgenomics.FgBioDef._
      val iterator = lines.bufferBetter

      val header = iterator.head
      if (header.startsWith("geneName")) {
        iterator
      }
      else if (header.startsWith("GENE_NAME")) {
        // convert the Picard-style header
        Iterator(header.split('\t').map(RefFlatSource.PicardHeaderMap(_)).mkString("\t")) ++ iterator.drop(1)
      }
      else {
        Iterator(RefFlatSource.Header.mkString("\t")) ++ iterator
      }
    }


    // Each line is a gene and transcript.  We want to group all transcripts from the same gene.
    var numDifferentChromosomes = 0
    var numDifferentStrands     = 0
    var numTranscripts          = 0
    val _genes = new DelimitedDataParser(lines=_lines, delimiter='\t').flatMap { row =>
      val geneName       = row[String]("geneName")
      val transcriptName = row[String]("name")
      val contig         = row[String]("chrom")
      val strand         = row[String]("strand")
      val exonCount      = row[Int]("exonCount")
      val exonStarts     = row[String]("exonStarts").split(',').filter(_.nonEmpty).map(_.toInt)
      val exonEnds       = row[String]("exonEnds").split(',').filter(_.nonEmpty).map(_.toInt)

      if (dict.exists(_.getSequence(contig) == null)) { // skip unrecognized sequences
        None
      }
      else {
        require(exonStarts.length == exonCount, s"Number of exonStarts does not equal exonCount for $geneName/$transcriptName")
        require(exonEnds.length == exonCount, s"Number of exonEnds does not equal exonCount for $geneName/$transcriptName")

        // Convert from 0-based half-open to 1-based inclusive intervals.
        val transcript = Transcript(
          name     = transcriptName,
          start    = row[Int]("txStart") + 1,
          end      = row[Int]("txEnd"),
          cdsStart = row[Int]("cdsStart") + 1,
          cdsEnd   = row[Int]("cdsEnd"),
          exons    = exonStarts.zip(exonEnds).map { case (s, e) => Exon(start=s + 1, end=e) }
        )

        val gene = Gene(
          contig         = contig,
          start          = -1,
          end            = -1,
          negativeStrand = strand == "-",
          name           = geneName,
          transcripts    = Seq(transcript)
        )

        Some(gene)
      }
    }
      .toSeq
      .groupBy(_.name)
      .map { case (name, _genes) =>
        val contig         = _genes.head.contig
        val negativeStrand = _genes.head.negativeStrand

        val transcripts = _genes.flatMap { gene =>
          require(gene.transcripts.length == 1, s"Found more than one transcript for gene ${gene.name}.")
          if (gene.contig != contig) {
            val transcript = gene.head
            numDifferentChromosomes += 1
            logger.info(s"Filtering out transcript '${transcript.name}' for gene ${gene.name}' due to being on a different chromosomes (${gene.contig}) than the first transcript seen ($contig).")
            None
          }
          else if (gene.negativeStrand != negativeStrand) {
            val transcript = gene.head
            val originalStrand = if (negativeStrand)      "-" else "+"
            val currentStrand  = if (gene.negativeStrand) "-" else "+"
            numDifferentStrands += 1
            logger.info(s"Filtering out transcript '${transcript.name}' for gene ${gene.name}' due to being on a different strand ($currentStrand) than the first transcript seen ($originalStrand).")
            None
          }
          else {
            numTranscripts += 1
            gene.transcripts
          }
        }

        Gene(
          contig         = contig,
          start          = transcripts.map(_.start).min,
          end            = transcripts.map(_.end).min,
          negativeStrand = negativeStrand,
          name           = name,
          transcripts    = transcripts
        )
      }

    logger.info(f"Filtered out $numDifferentChromosomes out of $numTranscripts (${numDifferentChromosomes/numTranscripts.toDouble*100}%.2f%%) transcript(s) due to being on a different chromosome.")
    logger.info(f"Filtered out $numDifferentStrands out of $numTranscripts (${numDifferentStrands/numTranscripts.toDouble*100}%.2f%%) transcript(s) due to being on a different strands.")

    _genes
  }

  def iterator: Iterator[Gene] = this.genes.toIterator

  /** Closes the underlying reader; only necessary if EOF hasn't been reached. */
  override def close(): Unit = this.source.foreach(_.close())
}
