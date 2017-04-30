/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.io.Io
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.{ProgressLogger, _}
import htsjdk.samtools.util._

@clp(description =
  """
    |Extracts unique molecular indexes from reads in a BAM file into tags.
    |
    |Currently only unmapped reads are supported.
    |
    |Only template bases will be retained as read bases (stored in the `SEQ` field) as specified by the read structure.
    |
    |A read structure should be provided for each read of a template.  For example, paired end reads should have two
    |read structures specified.  The tags to store the molecular indices will be associated with the molecular index
    |segment(s) in the read structure based on the order specified.  If only one molecular index tag is given, then the
    |molecular indices will be concatenated and stored in that tag. Otherwise the number of molecular indices in the
    |read structure should match the number of tags given. In the resulting BAM file each end of a pair will contain
    |the same molecular index tags and values. Additionally, when multiple molecular indices are present the
    |`--single-tag` option may be used to write all indices, concatenated, to a single tag in addition to the tags
    |specified in `--molecular-index-tags`.
    |
    |Optionally, the read names can be annotated with the molecular indices directly.  In this case, the read name
    |will be formatted `<NAME>+<UMIs1><UMIs2>` where `<UMIs1>` is the concatenation of read one's molecular indices.
    |Similarly for `<UMIs2>`.
    |
    |Mapping information will not be adjusted, as such, this tool should not be used on reads that have been mapped since
    |it will lead to an BAM with inconsistent records.
    |
    |The read structure describes the structure of a given read as one or more read segments. A read segment describes
    |a contiguous stretch of bases of the same type (ex. template bases) of some length and some offset from the start
    |of the read.  Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files.
    |Four kinds ofoperators are recognized:
    |
    |1. `T` identifies a template read
    |2. `B` identifies a sample barcode read
    |3. `M` identifies a unique molecular index read
    |4. `S` identifies a set of bases that should be skipped or ignored
    |
    |The last `<number><operator>` pair may be specified using a '+' sign instead of number to denote "all remaining
    |bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length.
    |
    |An example would be `10B3M7S100T` which describes 120 bases, with the first ten bases being a sample barcode,
    |bases 11-13 being a molecular index, bases 14-20 ignored, and bases 21-120 being template bases. See
    |[Read Structures](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) for more information.
  """,
  group = ClpGroups.SamOrBam)
class ExtractUmisFromBam
( @arg(flag='i', doc = "Input BAM file.")                                      val input: PathToBam,
  @arg(flag='o', doc = "Output BAM file.")                                     val output: PathToBam,
  @arg(flag='r', doc = "The read structure, one per read in a template.")      val readStructure: Seq[ReadStructure],
  @deprecated("Use molecular-index-tags instead.", since="0.1.3")
  @arg(flag='b', doc = "**[DEPRECATED]** SAM tags in which to store the molecular barcodes (one-per segment).",
    mutex=Array("molecularIndexTags"), minElements=0) val molecularBarcodeTags: Seq[String] = Seq.empty,
  @arg(flag='t', doc = "SAM tag(s) in which to store the molecular indices.", mutex=Array("molecularBarcodeTags"), minElements=0)
                                                                                 val molecularIndexTags: Seq[String] = Seq.empty,
  @arg(flag='s', doc = "Single tag into which to concatenate all molecular indices.") val singleTag: Option[String] = None,
  @arg(flag='a', doc = "Annotate the read names with the molecular indices. See usage for more details.") val annotateReadNames: Boolean = false,
  @arg(flag='c', doc = "The SAM tag with the position in read to clip adapters (e.g. `XT` as produced by Picard's `MarkIlluminaAdapters`).") val clippingAttribute: Option[String] = None
) extends FgBioTool with LazyLogging {

  val progress = ProgressLogger(logger, verb="written", unit=5e6.toInt)
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  val (rs1, rs2) = readStructure match {
    case Seq(r1, r2) => (r1.withVariableLastSegment, Some(r2.withVariableLastSegment))
    case Seq(r1)     => (r1.withVariableLastSegment, None)
    case Seq() => invalid("No read structures given")
    case _     => invalid("More than two read structures given")
  }

  // This can be removed once the @deprecated molecularBarcodeTags is removed
  private val perIndexTags = if (molecularIndexTags.nonEmpty) molecularIndexTags else molecularBarcodeTags
  if (perIndexTags.isEmpty) invalid("At least one molecular-index-tag must be specified.")

  // validate the read structure versus the molecular index tags
  {
    // create a read structure for the entire template
    val rsMolecularBarcodeSegmentCount = readStructure.map(_.molecularBarcodeSegments.size).sum
    // make sure each tag is of length 2
    perIndexTags.foreach(tag => if (tag.length != 2) invalid("SAM tags must be of length two: " + tag))
    // ensure we either have one tag, or we have the same # of tags as molecular indices in the read structure.
    if (perIndexTags.size > 1 && rsMolecularBarcodeSegmentCount != perIndexTags.size) {
      invalid("Either a single SAM tag, or the same # of SAM tags as molecular indices in the read structure,  must be given.")
    }
  }

  // split them into to tags for segments in read 1 vs read 2
  private val (molecularIndexTags1, molecularIndexTags2) = perIndexTags match {
    case Seq(tag) => (Seq[String](tag), Seq[String](tag))
    case tags =>
      val numMolecularIndicesRead1 = rs1.molecularBarcodeSegments.length
      (tags.slice(0, numMolecularIndicesRead1), tags.slice(numMolecularIndicesRead1, tags.length))
  }

  // Verify that if a single tag was specified that it is valid and not also contained in the per-index tags
  singleTag.foreach { tag =>
    if (tag.length != 2) invalid(s"All tags must be two characters: ${tag}")
    if (perIndexTags.contains(tag)) invalid(s"$tag specified as single-tag and also per-index tag")
  }

  override def execute(): Unit = {
    assert(molecularIndexTags1.length == rs1.molecularBarcodeSegments.length)
    rs2.foreach(rs => assert(molecularIndexTags2.length == rs.molecularBarcodeSegments.length))

    val in       = SamSource(input)
    val iterator = in.iterator.buffered
    val out      = SamWriter(output, in.header)

    while (iterator.hasNext) {
      val r1 = iterator.next()

      if (r1.paired) {
        // get the second read and make sure there's a read structure.
        if (!iterator.hasNext) fail(s"Could not find mate for read ${r1.name}")
        if (rs2.isEmpty) fail(s"Missing read structure for read two (required for paired end data).")
        val r2  = iterator.next()

        // verify everything is in order for paired end reads
        Seq(r1, r2).foreach(r => {
          if (!r.paired) fail(s"Read ${r.name} was not paired")
          if (r.mapped)  fail(s"Read ${r.name} was not unmapped")
        })
        if (!r1.firstOfPair)  fail(s"Read ${r1.name} was not the first end of a pair")
        if (!r2.secondOfPair) fail(s"Read ${r2.name} was not the second end of a pair")
        if (!r1.name.equals(r2.name)) fail(s"Read names did not match: '${r1.name}' and '${r2.name}'")

        // now do some work
        val bases1 = ExtractUmisFromBam.annotateRecord(record=r1, readStructure=rs1, molecularIndexTags=molecularIndexTags1, clippingAttribute=clippingAttribute)
        val bases2 = ExtractUmisFromBam.annotateRecord(record=r2, readStructure=rs2.get, molecularIndexTags=molecularIndexTags2, clippingAttribute=clippingAttribute)
        if (annotateReadNames) {
          // Developer Note: the delimiter must be an ascii character less than what is usually in the read names.  For
          // example "|" doesn't work.  I am not completely sure why.
          r1.name = r1.name + "+" + bases1 + bases2
          r2.name = r2.name + "+" + bases1 + bases2
        }
        assert(r1.name.equals(r2.name))

        // If we have duplicate tags, then concatenate them in the same order across the read pair.
        val tagAndValues = (molecularIndexTags1.map { tag => (tag, r1[String](tag)) } ++ molecularIndexTags2.map { tag => (tag, r2[String](tag)) }).toList
        tagAndValues.groupBy(_._1).foreach { case (tag, tuples) =>
          val attr = tuples.map(_._2).mkString(ExtractUmisFromBam.UmiDelimiter)
          r1(tag) = attr
          r2(tag) = attr
        }

        // If we have a single-tag, then also output values there
        singleTag.foreach { tag =>
          val value = tagAndValues.map { case (t,v) => v }.mkString(ExtractUmisFromBam.UmiDelimiter)
          r1(tag) = value
          r2(tag) = value
        }

        val recs = Seq(r1, r2)
        out ++= recs
        recs.foreach(progress.record)
      }
      else {
        // verify everything is in order for single end reads
        if (r1.mapped) fail(s"Read ${r1.name} was not unmapped")

        // now do some work
        val bases1 = ExtractUmisFromBam.annotateRecord(record=r1, readStructure=rs1, molecularIndexTags=molecularIndexTags1, clippingAttribute=clippingAttribute)
        if (annotateReadNames) {
          // Developer Note: the delimiter must be an ascii character less than what is usually in the read names.  For
          // example "|" doesn't work.  I am not completely sure why.
          r1.name = r1.name + "+" + bases1
        }

        // If we have duplicate tags, then concatenate them in the same order across the read
        val tagAndValues = molecularIndexTags1.map { tag => (tag, r1[String](tag)) }
        tagAndValues.groupBy(_._1).foreach { case (tag, tuples) =>
          val attr = tuples.map(_._2).mkString(ExtractUmisFromBam.UmiDelimiter)
          r1(tag) = attr
        }

        out += r1
        progress.record(r1)
      }
    }

    out.close()
    CloserUtil.close(iterator)
  }
}

object ExtractUmisFromBam {
  val UmiDelimiter = "-"

  /**
    * Extracts bases associated with molecular indices and adds them as SAM tags.  The read's bases will only contain
    * template bases as defined in the given read structure.
    */
  private[umi] def annotateRecord(record: SamRecord,
                                  readStructure: ReadStructure,
                                  molecularIndexTags: Seq[String],
                                  clippingAttribute: Option[String] = None): String = {
    val bases = record.basesString
    val qualities = record.qualsString
    // get the bases associated with each segment
    val readStructureBases = readStructure.extract(bases, qualities)
    // get the molecular index segments
    val molecularIndexBases = readStructureBases.filter(_.kind == SegmentType.MolecularBarcode).map(_.bases)

    // set the index tags
    // TODO: when we remove the deprecated molecularBarcodeTags option, consider whether or not we still
    //       need to have support for specifying a single tag via molecularIndexTags.
    molecularIndexTags match {
      case Seq(tag) => record(tag) = molecularIndexBases.mkString(UmiDelimiter)
      case _ =>
        if (molecularIndexTags.length < molecularIndexBases.length) throw new IllegalStateException("Found fewer molecular index SAM tags than molecular indices in the read structure.")
        else if (molecularIndexTags.length > molecularIndexBases.length) throw new IllegalStateException("Found fewer molecular indices in the read structure than molecular index SAM tags.")
        molecularIndexTags.zip(molecularIndexBases).foreach { case (tag, b) => record(tag) =  b }
    }

    // keep only template bases and qualities in the output read
    val basesAndQualities = readStructureBases.filter(_.kind == SegmentType.Template)

    // update any clipping information
    updateClippingInformation(record=record, clippingAttribute=clippingAttribute, readStructure=readStructure)
    record.bases = basesAndQualities.map(_.bases).mkString
    record.quals = basesAndQualities.map(_.quals).mkString

    // return the concatenation of the molecular indices in sequencing order
    molecularIndexBases.mkString
  }


  /**
    * Update the clipping information produced by Picard's MarkIlluminaAdapters to account for any non-template bases
    * that will be removed from the read.
    */
  private[umi] def updateClippingInformation(record: SamRecord,
                                             clippingAttribute: Option[String],
                                             readStructure: ReadStructure): Unit = {
    clippingAttribute.map(tag => (tag, record.get[Int](tag))) match {
      case None => Unit
      case Some((tag, None)) => Unit
      case Some((tag, Some(clippingPosition))) =>
        val newClippingPosition = readStructure.takeWhile(_.offset < clippingPosition).filter(_.kind == SegmentType.Template).map { t =>
          if (t.length.exists(l => t.offset + l < clippingPosition)) t.length.get
          else clippingPosition - t.offset
        }.sum
        record(tag) = newClippingPosition
    }
  }
}
