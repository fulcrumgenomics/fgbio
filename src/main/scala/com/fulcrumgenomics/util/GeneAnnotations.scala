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
                        exons: Seq[Exon],
                        features: Seq[Feature] = Seq(),
                        kind: Option[String] = None) extends Locatable {
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

    /** The order in which features appear in the transcripts */
    def transcriptOrderFeatures: Iterator[Feature] = features.iterator

    /** The order in which features appear in the genome */
    def genomicOrderFeatures: Iterator[Feature] = features.sortBy { feature => (feature.start, feature.end) }.iterator

    /** The features grouped by `kind` */
    def featuresByKind: Map[String, Seq[Feature]] = features.groupBy(f => f.kind)
  }

  /** Defines an exonic sequence within a transcript. */
  case class Exon(start: Int, end: Int) {
    require(start <= end, "start is greater than end when creating an exon")
    def length: Int = end - start + 1
  }

  /** Defines a generic genomic feature within a transcript. */
  case class Feature(id: String, start: Int, end: Int, kind: String)


  /** Trait that gene biotypes will extends */
  sealed trait GeneBiotype extends EnumEntry {
    /** The biotype name in the GFF */
    def key: String
    /** True if this biotype represents a protein coding feature */
    def isCoding: Boolean
  }

  /** Enum to represent gene biotypes in a GFF */
  object GeneBiotype extends FgBioEnum[GeneBiotype] {
    /** All values of the GeneBiotype Enum */
    def values: IndexedSeq[GeneBiotype] = findValues

    /** Allows GeneBiotypes to be built from the Enum name of the biotype name from the GFF */
    override def apply(str: String): GeneBiotype = values.find(_.key == str).getOrElse(super.apply(str))

    // antisense_RNA -> exon
    case object AntisenseRna  extends GeneBiotype { val key: String = "antisense_RNA";  val isCoding = false } // TODO: check this

    // C_gene_segment -> {exon, CDS}
    case object CRegion       extends GeneBiotype { val key: String = "C_region";       val isCoding = true }

    // D_gene_segment -> {exon, CDS}
    case object DSegment      extends GeneBiotype { val key: String = "D_segment";      val isCoding = true }

    // guide_RNA -> exon
    case object GuideRna      extends GeneBiotype { val key: String = "guide_RNA";      val isCoding = false }

    // J_gene_segment -> {exon, CDS}
    case object JSegment      extends GeneBiotype { val key: String = "J_segment";      val isCoding = true }

    // lnc_RNA -> exon
    case object LncRna        extends GeneBiotype { val key: String = "lncRNA";         val isCoding = false }

    // primary_transcript -> {exon, miRNA}  <-- odd
    case object MiRna         extends GeneBiotype { val key: String = "miRNA";          val isCoding = false }

    // transcript -> exon                   <-- odd
    case object MiscRna       extends GeneBiotype { val key: String = "misc_RNA";       val isCoding = false } // TODO: check this

    // *_feature or nothing                 <-- odd
    case object Other         extends GeneBiotype { val key: String = "other";          val isCoding = false } // TODO: check this

    // { transcript(for noncoding txs), mRNA(for coding txs) } -> exon         <-- odd
    case object ProteinCoding extends GeneBiotype { val key: String = "protein_coding"; val isCoding = true  }

    // RNase_MRP_RNA -> exon
    case object RNaseMrpRna   extends GeneBiotype { val key: String = "RNase_MRP_RNA";  val isCoding = false }

    // RNase_P_RNA -> exon
    case object RNasePRna     extends GeneBiotype { val key: String = "RNase_P_RNA";    val isCoding = false }

    // rRNA -> exon
    case object RRna          extends GeneBiotype { val key: String = "rRNA";           val isCoding = false }

    // scRNA -> exon
    case object ScRna         extends GeneBiotype { val key: String = "scRNA";          val isCoding = false }

    // snRNA -> exon
    case object SnRna         extends GeneBiotype { val key: String = "snRNA";          val isCoding = false }

    // snoRNA -> exon
    case object SnoRna        extends GeneBiotype { val key: String = "snoRNA";         val isCoding = false }

    // tRNA -> exon
    case object TRna          extends GeneBiotype { val key: String = "tRNA";           val isCoding = false }

    // telomerase_RNA -> exon
    case object TelomeraseRna extends GeneBiotype { val key: String = "telomerase_RNA"; val isCoding = false }

    // vault_RNA -> exon
    case object VaultRna      extends GeneBiotype { val key: String = "vault_RNA";      val isCoding = false }

    // V_segment -> exon
    case object VSegment      extends GeneBiotype { val key: String = "V_segment";      val isCoding = false }

    // Y_RNA -> exon
    case object YRna          extends GeneBiotype { val key: String = "Y_RNA";          val isCoding = false }
  }
}
