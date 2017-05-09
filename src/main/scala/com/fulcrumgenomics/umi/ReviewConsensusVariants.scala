/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.umi

import java.nio.file.Path
import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.BaseCounts
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.umi.ReviewConsensusVariants._
import com.fulcrumgenomics.util.{Io, Metric}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.IOUtil.{VCF_EXTENSIONS => VcfExtensions}
import htsjdk.samtools.util._
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.JavaConverters._
import scala.collection.mutable
import scala.collection.mutable.ListBuffer

object ReviewConsensusVariants {

  /**
    * Metric class that is used to output detailed information on each consensus read that carries a variant allele
    * @param chrom the chromosome on which the variant exists
    * @param pos the position of the variant
    * @param ref the reference allele at the position
    * @param genotype the genotype of the sample in question
    * @param filters the set of filters applied to the variant in the VCF
    * @param A the count of A observations at the variant locus
    * @param C the count of C observations at the variant locus
    * @param G the count of G observations at the variant locus
    * @param T the count of T observations at the variant locus
    * @param N the count of N observations at the variant locus
    * @param consensus_read the consensus read name for which the following fields contain values
    * @param consensus_insert a description of the insert that generated the consensus read
    * @param consensus_call the base call from the consensus read
    * @param consensus_qual the quality score from the consensus read
    * @param a the number of As in raw-reads contributing to the consensus base call at the variant site
    * @param c the number of Cs in raw-reads contributing to the consensus base call at the variant site
    * @param g the number of Gs in raw-reads contributing to the consensus base call at the variant site
    * @param t the number of Ts in raw-reads contributing to the consensus base call at the variant site
    * @param n the number of Ns in raw-reads contributing to the consensus base call at the variant site
    */
  case class ConsensusVariantReviewInfo
  ( // First set is all the same for each variant
    chrom: String, pos: Int, ref: String, genotype: String, filters: String, A: Int, C: Int, G: Int, T: Int, N: Int,
    // Second set is per-consensus reads
    consensus_read: String, consensus_insert: String, consensus_call: String, consensus_qual: Int, a: Int, c: Int, g: Int, t: Int, n: Int
  ) extends Metric

  /** Extracts the molecular identifier minus any trailing /A etc. */
  private[umi] def toMi(r: SAMRecord): String = Option(r.getStringAttribute(ConsensusTags.MolecularId)) match {
    case Some(mi) =>
      val slash = mi.lastIndexOf('/')
      if (slash > 0) mi.substring(0, slash) else mi
    case None     =>
      throw new IllegalStateException(s"${r.getReadName} did not have a value for tag ${ConsensusTags.MolecularId}")
  }

  /** Generate a /1 or /2 suffix based on whether a read is first or second of pair or unpaired. */
  private[umi] def readNumberSuffix(r : SAMRecord): String = if (r.getReadPairedFlag && r.getSecondOfPairFlag) "/2" else "/1"

  /**
    * Generates a string of the format 'chr1:500-650 | F1R2' for the insert represented by a read.
    * Only generates strings for FR pairs, all other reads get "NA".
    */
  private[umi] def toInsertString(r: SAMRecord): String = {
    if (!r.getReadPairedFlag || r.getReadUnmappedFlag || r.getMateUnmappedFlag ||
      r.getReferenceIndex != r.getMateReferenceIndex ||
      SamPairUtil.getPairOrientation(r) != PairOrientation.FR) {
      "NA"
    }
    else {
      val outer   = if (r.getReadNegativeStrandFlag) r.getAlignmentEnd else r.getAlignmentStart
      val isize   = r.getInferredInsertSize
      val Seq(start, end) = Seq(outer, outer + isize + (if (isize < 0) 1 else -1)).sorted

      val pairing = (r.getFirstOfPairFlag, start == outer) match {
        case (true, true)   => "F1R2"
        case (true, false)  => "F2R1"
        case (false, true)  => "F2R1"
        case (false, false) => "F1R2"
      }

      s"${r.getReferenceName}:${start}-${end} | ${pairing}"
    }
  }
}

@sopt.clp(group=ClpGroups.Umi, description=
  """
    |Extracts data to make reviewing of variant calls from consensus reads easier. Creates
    |a list of variant sites from the input VCF (SNPs only) or IntervalList then extracts all
    |the consensus reads that do not contain a reference allele at the the variant sites, and
    |all raw reads that contributed to those consensus reads.  This will include consensus
    |reads that carry the alternate allele, a third allele, a no-call or a spanning
    |deletion at the variant site.
    |
    |Reads are correlated between consensus and grouped BAMs using a molecule ID stored
    |in an optional attribute, MI by default.  In order to support paired molecule IDs
    |where two or more molecule IDs are related (e.g. see the Paired assignment strategy
    |in GroupReadsByUmi) the molecule ID is truncated at the last '/' if present
    |(e.g. "1/A" => "1" and "2" => "2").
    |
    |Both input BAMs must be coordinate sorted and indexed.
    |
    |A pair of output BAMs named <output>.consensus.bam and <output>.grouped.bam are created
    |with the relevant reads from each input BAM, and a review file <output>.txt is
    |created.  The review file contains details on each variant position along with detailed
    |information on each consensus read that supports the variant.  If the sample-name argument
    |is supplied and the input is VCF, genotype information for that sample will be retrieved.
    |If the sample-name isn't supplied and the VCF contains only a single sample then those
    |genotypes will be used.
  """)
class ReviewConsensusVariants
( @arg(flag='i', doc="Input VCF or IntervalList of variant locations.") val input : FilePath,
  @arg(flag='s', doc="Name of the sample being reviewed.") val sample: Option[String] = None,
  @arg(flag='c', doc="BAM file of consensus reads used to call variants.")  val consensusBam : PathToBam,
  @arg(flag='g', doc="BAM file of grouped raw reads used to build consensuses.") val groupedBam : PathToBam,
  @arg(flag='r', doc="Reference fasta file.") val ref: PathToFasta,
  @arg(flag='o', doc="Basename of output files to create.") val output : PathPrefix,
  @arg(flag='N', doc="Ignore N bases in the consensus reads.") val ignoreNsInConsensusReads: Boolean = false,
  @arg(flag='m', doc="Only output detailed information for variants at maf and below.") val maf: Double = 0.05
)extends FgBioTool with LazyLogging {
  Io.assertReadable(Seq(input, consensusBam, groupedBam, ref))
  Io.assertCanWriteFile(output)
  private val refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
  private val dict = refFile.getSequenceDictionary

  /** Simple case class to hold the relevant information from a variant context for easy access. */
  private[umi] case class Variant(chrom: String, start: Int, refBase: Char, genotype: Option[String], filters: Option[String]) {
    def toQueryInterval: QueryInterval = new QueryInterval(dict.getSequenceIndex(chrom), start, start)
    def toInterval: Interval = new Interval(chrom, start, start)
  }

  override def execute(): Unit = {
    def f(ext: String): Path = output.getParent.resolve(output.getFileName + ext)

    val consensusIn  = SamReaderFactory.make().open(consensusBam)
    val groupedIn    = SamReaderFactory.make().open(groupedBam)
    val factory      = new SAMFileWriterFactory().setCreateIndex(true)
    val consensusOut = factory.makeWriter(consensusIn.getFileHeader, true, f(".consensus.bam").toFile, ref.toFile)
    val groupedOut   = factory.makeWriter(groupedIn.getFileHeader,   true, f(".grouped.bam").toFile,   ref.toFile)

    // Pull out the consensus reads
    logger.info(s"Loading variants from ${input.toAbsolutePath}")
    val variants        = loadPositions(path=input, refFile=refFile)
    logger.info(s"Loaded ${variants.size} variant positions for review.")
    val variantsByChrom = variants.groupBy(v => v.chrom)
    val queries         = variants.map(_.toQueryInterval).toArray
    val sources         = mutable.HashSet[String]()

    logger.info(s"Extracting consensus reads from ${consensusBam.toAbsolutePath}")
    consensusIn.query(queries, false).filter(r => nonReferenceAtAnyVariant(r, variantsByChrom)).foreach { rec =>
      consensusOut.addAlignment(rec)
      sources.add(toMi(rec))
    }
    consensusIn.safelyClose()
    consensusOut.close()

    // And now the raw reads
    logger.info(s"Extracting raw/grouped reads from ${groupedBam.toAbsolutePath}")
    groupedIn.query(queries, false)
      .filter(rec => sources.contains(toMi(rec)))
      .foreach(groupedOut.addAlignment)
    groupedIn.safelyClose()
    groupedOut.close()

    // Generate the details output
    logger.info("Generating review file.")
    generateDetailsFile(consensusBam, f(".grouped.bam"), variants, f(".txt"))
  }

  /**
    * Generates the detailed consensus review file.
    *
    * @param consensusBam the original consensus bam (not the review consensus bam)
    * @param groupedBam the grouped bam (the review version)
    * @param variants the set of variants to report on
    * @param output the file to write
    */
  private def generateDetailsFile(consensusBam: PathToBam, groupedBam: PathToBam, variants: Seq[Variant], output: FilePath): Unit = {
    val consensusReader = SamReaderFactory.make().open(consensusBam)
    val groupedReader   = SamReaderFactory.make().open(groupedBam)
    val intervals = new IntervalList(consensusReader.getFileHeader)
    variants.foreach (v => intervals.add(v.toInterval))

    // Make locus iterators for both the consensus and grouped BAMs
    val Seq(consensusIterator, groupedIterator) = Seq(consensusReader, groupedReader).map { r =>
      val iterator = new SamLocusIterator(r, intervals, true)
      iterator.setEmitUncoveredLoci(true)
      iterator.setMappingQualityScoreCutoff(0)
      iterator.setQualityScoreCutoff(0)
      iterator.setSamFilters(Collections.emptyList())
      iterator.setIncludeIndels(true)
      iterator.iterator()
    }

    val writer = Metric.writer[ConsensusVariantReviewInfo](output)

    for (((variant, consensusPileup), groupedPileup) <- variants.iterator.zip(consensusIterator).zip(groupedIterator)) {
      require(variant.chrom == consensusPileup.getSequenceName, "Variants and consensus iterator out of sync.")
      require(variant.chrom == groupedPileup.getSequenceName,   "Variants and grouped iterator out of sync.")
      require(variant.start == consensusPileup.getPosition,     "Variants and consensus iterator out of sync.")
      require(variant.start == groupedPileup.getPosition,       "Variants and grouped iterator out of sync.")

      // Setup ready to iterate through the consensus reads
      val consensusCounts  = BaseCounts(consensusPileup.getRecordAndOffsets.asScala)
      val rawByMiAndReadNum = groupedPileup.getRecordAndOffsets.asScala.groupBy(r => toMi(r.getRecord) + readNumberSuffix(r.getRecord))

      consensusPileup.getRecordAndOffsets
        .filter(r => r.getReadBase.toChar.toUpper != variant.refBase)
        .filterNot(r => this.ignoreNsInConsensusReads && r.getReadBase.toChar.toUpper == 'N')
        .toSeq.sortBy(r => toMi(r.getRecord) + readNumberSuffix(r.getRecord))
        .foreach { c =>
          val mi = toMi(c.getRecord)
          val consensusReadName = c.getRecord.getReadName + readNumberSuffix(c.getRecord)
          val rawCounts = BaseCounts(rawByMiAndReadNum(mi + readNumberSuffix(c.getRecord)))

          val m =  ConsensusVariantReviewInfo(
            chrom = variant.chrom,
            pos   = variant.start,
            ref   = variant.refBase.toString,
            genotype = variant.genotype.getOrElse("NA"),
            filters  = variant.filters.getOrElse(""),
            A = consensusCounts.a,
            C = consensusCounts.c,
            G = consensusCounts.g,
            T = consensusCounts.t,
            N = consensusCounts.n,
            consensus_read = consensusReadName,
            consensus_insert = toInsertString(c.getRecord),
            consensus_call = c.getReadBase.toString.toUpperCase,
            consensus_qual = c.getBaseQuality,
            a = rawCounts.a,
            c = rawCounts.c,
            g = rawCounts.g,
            t = rawCounts.t,
            n = rawCounts.n
          )

          writer.write(m)
          writer.flush()
        }
    }

    writer.close()
  }

  /** Loads the set of locations for review. */
  private[umi] def loadPositions(path: FilePath, refFile: ReferenceSequenceFile): Seq[Variant] = {
    val buffer = ListBuffer[Variant]()

    if (VcfExtensions.exists(ext => path.getFileName.toString.endsWith(ext))) {
      val in = new VCFFileReader(input.toFile, false)
      val sample = this.sample.orElse {
        if (in.getFileHeader.getSampleNamesInOrder.size() == 1) Some(in.getFileHeader.getSampleNamesInOrder.get(0))
        else None
      }

      in.filter(_.isSNP).filter(v => mafFromGenotype(v, sample).forall(_ <= this.maf)).foreach { v => buffer += Variant(
        chrom    = v.getContig,
        start    = v.getStart,
        refBase  = v.getReference.getBases()(0).toChar.toUpper,
        genotype = sample.map(v.getGenotype(_).getGenotypeString),
        filters  = Some(v.getFilters.toSeq.sorted.mkString(","))
      )}
    }
    else {
      val list = IntervalList.fromFile(path.toFile).uniqued(false)
      for (i <- list; pos <- i.getStart to i.getEnd) {
        buffer += Variant(i.getContig, pos, refFile.getSubsequenceAt(i.getContig, pos, pos).getBases()(0).toChar.toUpper, None, None)
      }
    }

    buffer.toList
  }

  /** Attempts to determine the minor allele fraction from a VCF genotype object. */
  private def mafFromGenotype(ctx: VariantContext, sample: Option[String]): Option[Double] = sample match {
    case None    => None
    case Some(s) =>
      val gt = ctx.getGenotype(s)
      Option(gt.getExtendedAttribute("AF")) match {
        case Some(af) => Some(af.asInstanceOf[String].toDouble)
        case None =>
          val ad = gt.getAD
          if (ad == null || ad.isEmpty) None
          else Some (1 - ad(0) / ad.sum.toDouble)
      }
  }

  /**
    * Returns true if the read contains either a variant base, a no-call, or a deletion at one or more
    * of the variants that we are interested in.

    * @param rec the SAMRecord of interest
    * @param variantsByChromAndPos the set of variants we care about, grouped by chromosome
    * @return true if the read is non-reference at any variant base, otherwise false
    */
  private def nonReferenceAtAnyVariant(rec: SAMRecord, variantsByChromAndPos: Map[String, Seq[Variant]]): Boolean = {
    !rec.getReadUnmappedFlag && variantsByChromAndPos(rec.getContig).exists { v =>
      if (v.start >= rec.getAlignmentStart && v.start <= rec.getAlignmentEnd) {
        val readPos = rec.getReadPositionAtReferencePosition(v.start)
        if (readPos == 0) true
        else {
          val base = rec.getReadBases()(readPos - 1)
          !SequenceUtil.basesEqual(base, v.refBase.toByte) && (!this.ignoreNsInConsensusReads || !SequenceUtil.isNoCall(base))
        }
      }
      else {
        false
      }
    }
  }
}
