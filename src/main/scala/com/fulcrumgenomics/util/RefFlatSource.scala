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

import com.fulcrumgenomics.commons.CommonsDef.FilePath
import com.fulcrumgenomics.commons.util.{DelimitedDataParser, LazyLogging}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.util.GeneAnnotations.{Exon, Gene, GeneLocus, Transcript}

import scala.collection.mutable.ArrayBuffer
import scala.io.Source

object RefFlatSource {
  /** Creates a new refflat source from a sequence of lines. */
  def apply(lines: Seq[String], dict: Option[SequenceDictionary]): RefFlatSource = new RefFlatSource(lines.iterator, dict)

  /** Creates a new fasrefflattq source from an iterator of lines. */
  def apply(lines: Iterator[String], dict: Option[SequenceDictionary]): RefFlatSource = new RefFlatSource(lines, dict)

  /** Creates a new refflat source from an input stream. */
  def apply(stream: InputStream, dict: Option[SequenceDictionary]): RefFlatSource = new RefFlatSource(Source.fromInputStream(stream).getLines(), dict, Some(stream))

  /** Creates a new refflat source from a source. */
  def apply(source: Source, dict: Option[SequenceDictionary]): RefFlatSource = new RefFlatSource(source.getLines(), dict, Some(source))

  /** Creates a new refflat source from a File. */
  def apply(file: File, dict: Option[SequenceDictionary]): RefFlatSource = apply(path=file.toPath, dict=dict)

  /** Creates a new refflat source from a Path. */
  def apply(path: FilePath, dict: Option[SequenceDictionary] = None): RefFlatSource = apply(Io.toInputStream(path), dict)

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
                            dict: Option[SequenceDictionary] = None,
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
    val _genes = new DelimitedDataParser(lines=_lines, delimiter='\t').flatMap { row =>
      val geneName       = row.string("geneName")
      val transcriptName = row.string("name")
      val contig         = row.string("chrom").intern()
      val strand         = row.string("strand")
      val exonCount      = row.string("exonCount").toInt
      val exonStarts     = row.string("exonStarts").split(',').filter(_.nonEmpty).map(_.toInt)
      val exonEnds       = row.string("exonEnds").split(',').filter(_.nonEmpty).map(_.toInt)
      val isNegative     = strand == "-"

      if (dict.exists(!_.contains(contig))) { // skip unrecognized sequences
        None
      }
      else {
        require(exonStarts.length == exonCount, s"Number of exonStarts does not equal exonCount for $geneName/$transcriptName")
        require(exonEnds.length == exonCount, s"Number of exonEnds does not equal exonCount for $geneName/$transcriptName")

        // Convert from 0-based half-open to 1-based inclusive intervals.
        val exons = {
          val tmp = exonStarts.iterator.zip(exonEnds.iterator).map { case (s, e) => Exon(start=s + 1, end=e) }.toIndexedSeq
          if (isNegative) tmp.sortBy(e => -e.start) else tmp.sortBy(_.start)
        }

        // Detect where there is no coding region and set cds start and end appropriately
        val (cdsStart, cdsEnd) = (row[Int]("cdsStart"), row[Int]("cdsEnd")) match {
          case (s, e) if s == e => (None, None)
          case (s, e)           => (Some(s+1), Some(e))
        }

        val transcript = Transcript(
          name           = transcriptName,
          chrom          = contig,
          start          = row[Int]("txStart") + 1,
          end            = row[Int]("txEnd"),
          cdsStart       = cdsStart,
          cdsEnd         = cdsEnd,
          negativeStrand = isNegative,
          exons          = exons
        )

        val geneLocus = GeneLocus(Seq(transcript))
        Some(Gene(name=geneName, loci=Seq(geneLocus)))
      }
    }
      .toSeq
      .groupBy(_.name)
      .map { case (name, _genes) =>
        // Group the transcripts into groups that are:
        //   - On the same chromosome
        //   - One the same strand
        //   - Overlap at least one other transcript in the same group
        val groups = new ArrayBuffer[ArrayBuffer[Transcript]]()

        _genes.flatMap(_.iterator).sortBy(tx => -tx.lengthOnGenome).foreach { tx =>
          groups.find { group => group.exists(t => t.negativeStrand == tx.negativeStrand && t.overlaps(tx)) } match {
            case Some(group) => group += tx
            case None        => groups += ArrayBuffer(tx)
          }
        }

        Gene(name=name, loci=groups.map(g => GeneLocus(g.toSeq)).toSeq)
      }

    _genes
  }

  def iterator: Iterator[Gene] = this.genes.iterator

  /** Closes the underlying reader; only necessary if EOF hasn't been reached. */
  override def close(): Unit = this.source.foreach(_.close())
}
