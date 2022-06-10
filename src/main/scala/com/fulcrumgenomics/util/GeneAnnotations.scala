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

import htsjdk.samtools.util.{CoordMath, Locatable}
import com.fulcrumgenomics.FgBioDef._
import enumeratum.EnumEntry

/** Stores classes useful for storing annotation information for genes and their transcripts and exons. */
object GeneAnnotations {
  /**
    * A gene with a collection of gene loci.
    */
  case class Gene(name: String, loci: Seq[GeneLocus], biotype: Option[GeneBiotype] = None) extends Iterable[Transcript] {
    /** Iterate over all transcripts from all GeneLoci (each a collection of transcripts). */
    def iterator: Iterator[Transcript] = this.loci.iterator.flatten

    /** Returns the locus with the maximum transcript mappings. */
    def primaryLocus: GeneLocus = loci.maxBy(_.size)
  }


  /**
    * A gene locus representing a set of transcripts that are mapped to the same location on the genome.
    *
    * @param transcripts the set of transcripts at the locus
    */
  case class GeneLocus(transcripts: Seq[Transcript]) extends Locatable with Iterable[Transcript] {
    require(transcripts.nonEmpty, "GeneLocus cannot be generated from an empty set of transcripts.")
    require(transcripts.forall(t => t.chrom == transcripts.head.chrom))
    require(transcripts.forall(t => t.negativeStrand == transcripts.head.negativeStrand))

    val chrom: String           = transcripts.head.chrom
    val start: Int              = transcripts.minBy(_.start).start
    val end:   Int              = transcripts.maxBy(_.end).end
    val negativeStrand: Boolean = transcripts.head.negativeStrand
    def positiveStrand: Boolean = !negativeStrand

    // Locatable implementation
    override def getContig: String = chrom
    override def getStart: Int     = start
    override def getEnd: Int       = end

    override def iterator: Iterator[Transcript] = this.transcripts.iterator
  }


  /**
    * A transcript [mapping] associated with a given gene locus.  Contains zero or more exons.
    * The exons should be given in the order of transcription.
    */
  case class Transcript(name: String,
                        chrom: String,
                        start: Int,
                        end: Int,
                        cdsStart: Option[Int] = None,
                        cdsEnd: Option[Int]   = None,
                        negativeStrand: Boolean,
                        exons: Seq[Exon]) extends Locatable {
    if (exons.length > 1) {
      val exonsOverlap = genomicOrder.sliding(2).map { e => (e.head, e.last) }.exists { case (e1, e2) => CoordMath.overlaps(e1.start, e1.end, e2.start, e2.end) }
      require(!exonsOverlap, s"exons overlap for transcript: $name")
    }

    require(cdsStart.isDefined == cdsEnd.isDefined, s"cdsStart and cdsEnd must either both or neither be defined on tx: $name.")

    // Locatable implementation
    override def getContig: String = chrom
    override def getStart: Int = start
    override def getEnd: Int = end

    def positiveStrand: Boolean = !negativeStrand

    /** The order in which exons appear in the transcripts */
    def transcriptOrder: Iterator[Exon] = exons.iterator

    /** The order in which exons appear in the genome */
    def genomicOrder: Iterator[Exon] = exons.sortBy { exon => (exon.start, exon.end) }.iterator

    /** The length of the transcript mapping on the genome (including introns). */
    def lengthOnGenome: Int = end - start + 1

    /** The length of the transcribed product after introns are removed. */
    def transcribedLength: Int = exons.sumBy(_.length)

    /** True if the transcript is coding, false otherwise. */
    def isCoding: Boolean = cdsStart.isDefined
  }

  /** Defines an exonic sequence within a transcript. */
  case class Exon(start: Int, end: Int) {
    require(start <= end, "start is greater than end when creating an exon")
    def length: Int = end - start + 1
  }

  /** Trait that gene biotypes will extends */
  sealed trait GeneBiotype extends EnumEntry {
    /** The biotype name in the RefSeq GFF */
    def refSeqValue: String
    /** True if this biotype represents a protein coding feature */
    def isCoding: Boolean
  }

  /** Enum to represent gene biotypes.  The current list contains all values found in NCBI RefSeq GFFs circa 2021.
    * In future more and/or different values from other sources may be added. */
  object GeneBiotype extends FgBioEnum[GeneBiotype] {
    /** All values of the GeneBiotype Enum */
    def values: IndexedSeq[GeneBiotype] = findValues

    /** Allows GeneBiotypes to be built from the Enum name of the biotype name from the GFF */
    override def apply(str: String): GeneBiotype = values.find(_.refSeqValue == str).getOrElse(super.apply(str))

    case object AntisenseRna  extends GeneBiotype { val refSeqValue = "antisense_RNA";  val isCoding = false }
    case object CRegion       extends GeneBiotype { val refSeqValue = "C_region";       val isCoding = true }
    case object DSegment      extends GeneBiotype { val refSeqValue = "D_segment";      val isCoding = true }
    case object GuideRna      extends GeneBiotype { val refSeqValue = "guide_RNA";      val isCoding = false }
    case object JSegment      extends GeneBiotype { val refSeqValue = "J_segment";      val isCoding = true }
    case object LncRna        extends GeneBiotype { val refSeqValue = "lncRNA";         val isCoding = false }
    case object MiRna         extends GeneBiotype { val refSeqValue = "miRNA";          val isCoding = false }
    case object MiscRna       extends GeneBiotype { val refSeqValue = "misc_RNA";       val isCoding = false }
    case object NcRna         extends GeneBiotype { val refSeqValue = "ncRNA";          val isCoding = false }
    case object Other         extends GeneBiotype { val refSeqValue = "other";          val isCoding = false }
    case object ProteinCoding extends GeneBiotype { val refSeqValue = "protein_coding"; val isCoding = true  }
    case object RNaseMrpRna   extends GeneBiotype { val refSeqValue = "RNase_MRP_RNA";  val isCoding = false }
    case object RNasePRna     extends GeneBiotype { val refSeqValue = "RNase_P_RNA";    val isCoding = false }
    case object RRna          extends GeneBiotype { val refSeqValue = "rRNA";           val isCoding = false }
    case object ScRna         extends GeneBiotype { val refSeqValue = "scRNA";          val isCoding = false }
    case object SnRna         extends GeneBiotype { val refSeqValue = "snRNA";          val isCoding = false }
    case object SnoRna        extends GeneBiotype { val refSeqValue = "snoRNA";         val isCoding = false }
    case object SrpRna        extends GeneBiotype { val refSeqValue = "SRP_RNA";        val isCoding = false }
    case object TRna          extends GeneBiotype { val refSeqValue = "tRNA";           val isCoding = false }
    case object TelomeraseRna extends GeneBiotype { val refSeqValue = "telomerase_RNA"; val isCoding = false }
    case object VaultRna      extends GeneBiotype { val refSeqValue = "vault_RNA";      val isCoding = false }
    case object VSegment      extends GeneBiotype { val refSeqValue = "V_segment";      val isCoding = false }
    case object YRna          extends GeneBiotype { val refSeqValue = "Y_RNA";          val isCoding = false }
  }
}
