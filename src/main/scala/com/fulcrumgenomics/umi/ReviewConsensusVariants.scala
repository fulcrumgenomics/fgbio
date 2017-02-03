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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.Io
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.reference.{ReferenceSequenceFile, ReferenceSequenceFileFactory}
import htsjdk.samtools.util.IOUtil.{VCF_EXTENSIONS => VcfExtensions}
import htsjdk.samtools.util.{IntervalList, SequenceUtil}
import htsjdk.samtools.{QueryInterval, SAMFileWriterFactory, SAMRecord, SamReaderFactory}
import htsjdk.variant.vcf.VCFFileReader

import scala.collection.mutable
import scala.collection.mutable.ListBuffer

@clp(group=ClpGroups.Umi, description=
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
    |in GroupReadsByUmi) the molecule ID is truncated at the first '/' if present
    |(e.g. "1/A" => "1" and "2" => "2").
    |
    |Both input BAMs must be coordinate sorted and indexed.
    |
    |A pair of output BAMs named <output>.consensus.bam and <output>.grouped.bam
    |are created with the relevant reads from each input BAM.
  """)
class ReviewConsensusVariants
( @arg(flag="i", doc="Input VCF or IntervalList of variant locations.") val input : FilePath,
  @arg(flag="c", doc="BAM file of consensus reads used to call variants.")  val consensusBam : PathToBam,
  @arg(flag="g", doc="BAM file of grouped raw reads used to build consensuses.") val groupedBam : PathToBam,
  @arg(flag="r", doc="Reference fasta file.") val ref: PathToFasta,
  @arg(flag="o", doc="Basename of output files to create.") val output : PathPrefix,
  @arg(flag="t", doc="The SAM/BAM tag used to uniquely identify source molecules.") val tag: String = "MI",
  @arg(flag="N", doc="Ignore N bases in the consensus reads.") val ignoreNsInConsensusReads: Boolean = false
)extends FgBioTool with LazyLogging {
  Io.assertReadable(Seq(input, consensusBam, groupedBam, ref))
  Io.assertCanWriteFile(output)
  private val refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(ref)
  private val dict = refFile.getSequenceDictionary

  private[umi] case class LocAndRef(chrom: String, start: Int, refBase: Char) {
    def toQueryInterval: QueryInterval = new QueryInterval(dict.getSequenceIndex(chrom), start, start)
  }

  override def execute(): Unit = {
    def f(ext: String): Path = output.getParent.resolve(output.getFileName + ext)
    val consensusIn  = SamReaderFactory.make().open(consensusBam)
    val groupedIn    = SamReaderFactory.make().open(groupedBam)
    val factory      = new SAMFileWriterFactory().setCreateIndex(true)
    val consensusOut = factory.makeWriter(consensusIn.getFileHeader, true, f(".consensus.bam").toFile, ref.toFile)
    val groupedOut   = factory.makeWriter(groupedIn.getFileHeader, true, f(".grouped.bam").toFile, ref.toFile)

    // Pull out the consensus reads
    logger.info(s"Loading variants from ${input.toAbsolutePath}")
    val variants        = loadPositions(path=input, refFile=refFile)
    val variantsByChrom = variants.groupBy(v => v.chrom)
    val queries         = variants.map(_.toQueryInterval).toArray
    val sources         = mutable.HashSet[String]()

    logger.info(s"Extracting consensus reads from ${consensusBam.toAbsolutePath}")
    consensusIn.query(queries, false).filter(r => nonReferenceAtAnyVariant(r, variantsByChrom)).foreach { rec =>
      consensusOut.addAlignment(rec)
      Option(rec.getStringAttribute(tag)).map(mi => mi.takeWhile(_ != '/')).foreach(sources.add)
    }
    consensusIn.safelyClose()
    consensusOut.close()

    // And now the raw reads
    logger.info(s"Extracting raw/grouped reads from ${groupedBam.toAbsolutePath}")
    groupedIn.query(queries, false)
      .filter(rec => Option(rec.getStringAttribute(tag)).map(mi => mi.takeWhile(_ != '/')).exists(sources.contains))
      .foreach(groupedOut.addAlignment)
    groupedIn.safelyClose()
    groupedOut.close()
  }

  /** Loads the set of locations for review. */
  private[umi] def loadPositions(path: FilePath, refFile: ReferenceSequenceFile): Seq[LocAndRef] = {
    val buffer = ListBuffer[LocAndRef]()

    if (VcfExtensions.exists(ext => path.getFileName.toString.endsWith(ext))) {
      val in = new VCFFileReader(input.toFile, false)
      in.filter(_.isSNP).foreach(v => buffer += LocAndRef(v.getContig, v.getStart, v.getReference.getBases()(0).toChar.toUpper))
    }
    else {
      val list = IntervalList.fromFile(path.toFile).uniqued(false)
      for (i <- list; pos <- i.getStart to i.getEnd) {
        buffer += LocAndRef(i.getContig, pos, refFile.getSubsequenceAt(i.getContig, pos, pos).getBases()(0).toChar.toUpper)
      }
    }

    buffer.toList
  }

  /**
    * Returns true if the read contains either a variant base, a no-call, or a deletion at one or more
    * of the variants that we are interested in.

    * @param rec the SAMRecord of interest
    * @param variantsByChromAndPos the set of variants we care about, grouped by chromosome
    * @return true if the read is non-reference at any variant base, otherwise false
    */
  private def nonReferenceAtAnyVariant(rec: SAMRecord, variantsByChromAndPos: Map[String, Seq[LocAndRef]]): Boolean = {
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
