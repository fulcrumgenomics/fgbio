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
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.async.AsyncIterator
import com.fulcrumgenomics.commons.collection.BetterBufferedIterator
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools.{SAMFileHeader, SAMSequenceDictionary}
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor

/** Companion object that holds implementation methods for the ZipperBams command line tool. */
private[bam] object ZipperBams extends LazyLogging {
  /** Named sets of tags for reversing. */
  val TagSetsForReversing = Map(
    "Consensus" -> ConsensusTags.PerBase.TagsToReverse
  )

  /** Named sets of tags for reverse complementing. */
  val TagSetsForRevcomping = Map(
    "Consensus" -> ConsensusTags.PerBase.TagsToReverseComplement
  )

  object TagInfo {
    def apply(remove: IndexedSeq[String], reverse: IndexedSeq[String], revcomp: IndexedSeq[String]): TagInfo = {
      val r  = reverse.flatMap(t => TagSetsForReversing.getOrElse(t, Some(t)))
      val rc = revcomp.flatMap(t => TagSetsForRevcomping.getOrElse(t, Some(t)))
      new TagInfo(remove=remove, reverse=r, revcomp=rc)
    }
  }

  /** Class to hold info about all the extended tags/attrs we want to manipulate. */
  case class TagInfo private (remove: IndexedSeq[String], reverse: IndexedSeq[String], revcomp: IndexedSeq[String]) {
    val setToReverse: java.util.Set[String] = reverse.iterator.toJavaSet
    val setToRevcomp: java.util.Set[String] = revcomp.iterator.toJavaSet
    val hasRevsOrRevcomps: Boolean = reverse.nonEmpty || revcomp.nonEmpty
  }

  /** Checks the source's header declares it is queryname sorted or grouped and logs a warning if not. */
  def checkSort(source: SamSource, path: FilePath, name: String): Unit = {
    if (source.header.getSortOrder != SortOrder.queryname && source.header.getGroupOrder != GroupOrder.query) {
      logger.warning(s"${name} file ${path} does not appear to be queryname sorted or grouped.")
      logger.warning(s"Continuing, but your output may be incorrect.")
    }
  }

  /** Builds the header for the output BAM from the headers of the unmapped and mapped BAM.
    *
    * @param unmapped the header from the unmapped BAM
    * @param mapped the header from the mapped BAM
    * @param sort optionally, the sort order specified for the output BAM
    * @param dict the sequence dictionary of the reference genome to insert into the output header
    */
  def buildOutputHeader(unmapped: SAMFileHeader,
                        mapped: SAMFileHeader,
                        sort: Option[SamOrder],
                        dict: SAMSequenceDictionary): SAMFileHeader = {
    // Build a new header using the given sequence dictionary
    val header = new SAMFileHeader(dict)

    // Copy over comments, RGs and PGs, letting mapped override unmapped if there are conflicts
    Iterator(unmapped, mapped).foreach { old =>
      old.getComments.iterator().foreach(header.addComment)
      old.getReadGroups.iterator().foreach(header.addReadGroup)
      old.getProgramRecords.iterator().foreach(header.addProgramRecord)
    }

    // Set the sort and group order
    header.setSortOrder(unmapped.getSortOrder)
    header.setGroupOrder(unmapped.getGroupOrder)
    sort.foreach(s => s.applyTo(header))

    header
  }

  /** Constructs a template iterator that assumes the input is sorted or grouped, and turns it into an async iter. */
  def templateIterator(source: SamSource, bufferSize: Int): BetterBufferedIterator[Template] = {
    val recs      = source.iterator.bufferBetter
    val templates = new Iterator[Template] {
      override def hasNext: Boolean = recs.hasNext
      override def next(): Template = {
        try {
          Template(recs)
        }
        catch {
          case iae: IllegalArgumentException => throw new IllegalStateException(
            s"Error reading from ${source.toSamReader.getResourceDescription}: ${iae.getMessage}"
          )
        }
      }
    }

    if (bufferSize > 0) {
      AsyncIterator(templates, Some(bufferSize)).bufferBetter
    }
    else {
      templates.bufferBetter
    }
  }

  /** Merge information from the unmapped template into the mapped template. */
  def merge(unmapped: Template, mapped: Template, tagInfo: TagInfo): Template = {
    // Fix mate info first
    mapped.fixMateInfo()

    // Remove tags from the mapped records
    for (rec <- mapped.allReads; tag <- tagInfo.remove) rec.remove(tag)

    // Copy tags over from the unmapped record
    for (u <- unmapped.primaryReads) {
      // Get all the mapped records for the unmapped record
      val ms = if (u.unpaired || u.firstOfPair) mapped.allR1s.toSeq else mapped.allR2s.toSeq
      val hasPosReads = ms.exists(_.positiveStrand)
      val hasNegReads = ms.exists(_.negativeStrand)

      // Copy over pass/fail qc flag
      ms.foreach(_.pf = u.pf)

      // Copy over tags on any positive strand reads first
      if (hasPosReads) ms.iterator.filter(_.positiveStrand).foreach(m => copyTags(u, m))

      // Then reverse complement the unmapped read if we need to and do the negative strand reads
      if (hasNegReads) {
        if (tagInfo.hasRevsOrRevcomps) u.asSam.reverseComplement(tagInfo.setToRevcomp, tagInfo.setToReverse, !hasPosReads)
        ms.iterator.filter(_.negativeStrand).foreach(m => copyTags(u, m))
      }
    }

    mapped
  }

  /** Copies all tags from read to another, with the exception that PG is only copied if not already present on dest. */
  def copyTags(src: SamRecord, dest: SamRecord): Unit = {
    src.attributes.foreach { case (tag, value) =>
      if (tag != "PG" || !dest.contains("PG")) {
        dest(tag) = value
      }
    }
  }
}

@clp(group=ClpGroups.SamOrBam, description=
  """
    |Zips together an unmapped and mapped BAM to transfer metadata into the output BAM.
    |
    |Both the unmapped and mapped BAMs _must_ be a) queryname sorted or grouped (i.e. all records with the same
    |name are grouped together in the file), and b) have the same ordering of querynames.  If either of these are
    |violated the output is undefined!
    |
    |All tags present on the unmapped reads are transferred to the mapped reads.  The options `--tags-to-reverse`
    |and `--tags-to-revcomp` will cause tags on the unmapped reads to be reversed or reverse complemented before
    |being copied to reads mapped to the negative strand.  These options can take a mixture of two-letter tag names
    |and the names of tag sets, which will be expanded into sets of tag names.  Currently the only named tag set
    |is "Consensus" which contains all the per-base consensus tags produced by fgbio consensus callers.
    |
    |By default the mapped BAM is read from standard input (stdin) and the output BAM is written to standard
    |output (stdout). This can be changed using the `--input/-i` and `--output/-o` options.
    |
    |By default the output BAM file is emitted in the same order as the input BAMs.  This can be overridden
    |using the `--sort` option, though in practice it may be faster to do the following:
    |
    |```
    |fgbio --compression 0 ZipperBams -i mapped.bam -u unmapped.bam -r ref.fa | samtools sort -@ $(nproc)
    |```
  """)
class ZipperBams
( @arg(flag='i', doc="Mapped SAM or BAM.") val input: PathToBam = Io.StdIn,
  @arg(flag='u', doc="Unmapped SAM or BAM.") val unmapped: PathToBam,
  @arg(flag='r', doc="Path to the reference used in alignment. Must have accompanying .dict file.") val ref: PathToFasta,
  @arg(flag='o', doc="Output SAM or BAM file.") val output: PathToBam = Io.StdOut,
  @arg(doc="Tags to remove from the mapped BAM records.", minElements=0)
  val tagsToRemove: IndexedSeq[String] = IndexedSeq.empty,
  @arg(doc="Set of optional tags to reverse on reads mapped to the negative strand.", minElements=0)
  val tagsToReverse: IndexedSeq[String] = IndexedSeq.empty,
  @arg(doc="Set of optional tags to reverse complement on reads mapped to the negative strand.", minElements=0)
  val tagsToRevcomp: IndexedSeq[String] = IndexedSeq.empty,
  @arg(flag='s', doc="Sort the output BAM into the given order.") val sort: Option[SamOrder] = None,
  @arg(flag='b', doc="Buffer this many read-pairs while reading the input BAMs.") val buffer: Int = 5000
) extends FgBioTool {
  import ZipperBams._

  Io.assertReadable(ref)
  Io.assertReadable(input)
  Io.assertReadable(unmapped)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val unmappedSource = SamSource(unmapped)
    val mappedSource   = SamSource(input)
    checkSort(unmappedSource, unmapped, "unmapped")
    checkSort(mappedSource,  input, "input")

    val dict    = SAMSequenceDictionaryExtractor.extractDictionary(this.ref)
    val header  = buildOutputHeader(unmappedSource.header, mappedSource.header, this.sort, dict)
    val out     = SamWriter(output, header, sort=this.sort)
    val tagInfo = TagInfo(remove=tagsToRemove, reverse=tagsToReverse, revcomp=tagsToRevcomp)

    if (tagInfo.remove.nonEmpty)  logger.info(s"Tags for removal from mapped BAM: ${tagInfo.remove.mkString(", ")}")
    if (tagInfo.reverse.nonEmpty) logger.info(s"Tags being reversed: ${tagInfo.reverse.mkString(", ")}")
    if (tagInfo.revcomp.nonEmpty) logger.info(s"Tags being reverse complemented: ${tagInfo.revcomp.mkString(", ")}")

    val unmappedIter = templateIterator(unmappedSource, bufferSize=buffer)
    val mappedIter   = templateIterator(mappedSource, bufferSize=buffer)

    unmappedIter.foreach { template =>
      if (mappedIter.hasNext && mappedIter.head.name == template.name) {
        out ++= merge(unmapped=template, mapped=mappedIter.next(), tagInfo=tagInfo).allReads
      }
      else {
        logger.debug(s"Found an unmapped read with no corresponding mapped read: ${template.name}")
        out ++= template.allReads
      }
    }

    out.close()
    unmappedSource.safelyClose()
    mappedSource.safelyClose()
  }
}
