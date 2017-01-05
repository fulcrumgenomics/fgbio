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

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.ReadStructure.SubRead
import com.fulcrumgenomics.util.{ProgressLogger, _}
import com.fulcrumgenomics.FgBioDef.unreachable
import dagr.commons.CommonsDef.PathToBam
import dagr.commons.io.Io
import dagr.commons.util.LazyLogging
import dagr.sopt._
import dagr.sopt.cmdline.ValidationException
import htsjdk.samtools._
import htsjdk.samtools.util._

import scala.collection.JavaConversions._

@clp(description =
  """
    |Extracts unique molecular indexes from reads in a BAM file into tags.
    |
    |Currently only unmapped reads are supported.
    |
    |Only template bases will be retained as read bases (stored in the SEQ field) as specified by the read structure.
    |
    |A read structure should be provided for each read of a template.  For example, paired end reads should have two
    |read structures specified.  The tags to store the molecular indices will be associated with the molecular index
    |segment(s) in the read structure based on the order specified.  If only one molecular index tag is given, then the
    |molecular indices will be concatenated and stored in that tag. Otherwise the number of molecular indices in the
    |read structure should match the number of tags given. In the resulting BAM file each end of a pair will contain
    |the same molecular index tags and values. Additionally, when multiple molecular indices are present the
    |'single-tag' option may be used to write all indices, concatenated, to a single tag in addition to the tags
    |specified in 'molecular-index-tags'.
    |
    |Optionally, the read names can be annotated with the molecular indices directly.  In this case, the read name
    |will be formatted "<NAME>+<UMIs1><UMIs2>" where "<UMIs1>" is the concatenation of read one's molecular indices.
    |Similarly for "<UMIs2>".
    |
    |Mapping information will not be adjusted, as such, this tool should not be used on reads that have been mapped since
    |it will lead to an BAM with inconsistent records.
    |
    |The read structure describes the structure of a given read as one or more read segments. A read segment describes
    |a contiguous stretch of bases of the same type (ex. template bases) of some length and some offset from the start
    |of the read.  The following segment types are supported:
    |  - T: template bases
    |  - B: sample barcode bases
    |  - M: molecular index bases
    |  - S: bases to ignore
    |An example would be "10B3M7S100T" which describes 120 bases, with the first ten bases being a sample barcode,
    |bases 11-13 being a molecular index, bases 14-20 ignored, and bases 21-120 being template bases.
  """,
  group = ClpGroups.SamOrBam)
class ExtractUmisFromBam
( @arg(flag = "i", doc = "Input BAM file.")                                      val input: PathToBam,
  @arg(flag = "o", doc = "Output BAM file.")                                     val output: PathToBam,
  @arg(flag = "r", doc = "The read structure, one per read in a template.")      val readStructure: Seq[String],
  @deprecated @arg(flag = "b", doc = "[DEPRECATED] SAM tags in which to store the molecular barcodes (one-per segment).",
    mutex=Array("molecularIndexTags")) val molecularBarcodeTags: Seq[String] = Seq.empty,
  @arg(flag = "t", doc = "SAM tag(s) in which to store the molecular indices.", mutex=Array("molecularBarcodeTags"))
                                                                                 val molecularIndexTags: Seq[String] = Seq.empty,
  @arg(flag = "s", doc = "Single tag into which to concatenate all molecular indices.") val singleTag: Option[String] = None,
  @arg(flag = "a", doc = "Annotate the read names with the molecular indices. See usage for more details.") val annotateReadNames: Boolean = false,
  @arg(flag = "c", doc = "The SAM tag with the position in read to clip adapters (e.g. XT as produced by Picard's MarkIlluminaAdapters).") val clippingAttribute: Option[String] = None
) extends FgBioTool with LazyLogging {

  val progress = new ProgressLogger(logger, verb="written", unit=5e6.toInt)
  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  val (rs1, rs2) = readStructure match {
    case Seq(readStructure1, readStructure2) => (ReadStructure(readStructure1), Some(ReadStructure(readStructure2)))
    case Seq(readStructure1) => (ReadStructure(readStructure1), None)
    case Seq() => invalid("No read structures given")
    case _     => invalid("More than two read structures given")
  }

  // This can be removed once the @deprecated molecularBarcodeTags is removed
  private val perIndexTags = if (molecularIndexTags.nonEmpty) molecularIndexTags else molecularBarcodeTags

  // validate the read structure versus the molecular barcode tags
  {
    // create a read structure for the entire template
    val rs = ReadStructure(readStructure.mkString(""))
    // make sure each tag is of length 2
    perIndexTags.foreach(tag => if (tag.length != 2) invalid("SAM tags must be of length two: " + tag))
    // ensure we either have one tag, or we have the same # of tags as molecular barcodes in the read structure.
    if (perIndexTags.size > 1 && rs.molecularBarcode.size != perIndexTags.size) {
      invalid("Either a single SAM tag, or the same # of SAM tags as molecular barcodes in the read structure,  must be given.")
    }
  }

  // split them into to tags for segments in read 1 vs read 2
  private val (molecularIndexTags1, molecularIndexTags2) = perIndexTags match {
    case Seq(tag) => (Seq[String](tag), Seq[String](tag))
    case tags =>
      val numMolecularBarcodesRead1 = rs1.molecularBarcode.length
      (tags.subList(0, numMolecularBarcodesRead1).toSeq, tags.subList(numMolecularBarcodesRead1, tags.length).toSeq)
  }

  // Verify that if a single tag was specified that it is valid and not also contained in the per-index tags
  singleTag.foreach { tag =>
    if (tag.length != 2) invalid(s"All tags must be two characters: ${tag}")
    if (perIndexTags.contains(tag)) invalid(s"$tag specified as single-tag and also per-index tag")
  }

  override def execute(): Unit = {
    assert(molecularIndexTags1.length == rs1.count {
      case m: MolecularBarcode => true
      case _ => false
    })
    rs2.foreach { rs =>
      assert(molecularIndexTags2.length == rs.count {
        case m: MolecularBarcode => true
        case _ => false
      })
    }
    val in: SamReader = SamReaderFactory.make.open(input.toFile)
    val iterator = in.iterator().buffered
    val out: SAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(in.getFileHeader, true, output.toFile)

    while (iterator.hasNext) {
      val r1 = iterator.next()

      if (r1.getReadPairedFlag) {
        // get the second read and make sure there's a read structure.
        if (!iterator.hasNext) fail(s"Could not find mate for read ${r1.getReadName}")
        if (rs2.isEmpty) fail(s"Missing read structure for read two (required for paired end data).")
        val r2  = iterator.next()

        // verify everything is in order for paired end reads
        Seq(r1, r2).foreach(r => {
          if (!r.getReadPairedFlag) fail(s"Read ${r.getReadName} was not paired")
          if (!r.getReadUnmappedFlag) fail(s"Read ${r.getReadName} was not unmapped")
        })
        if (!r1.getFirstOfPairFlag) fail(s"Read ${r1.getReadName} was not the first end of a pair")
        if (!r2.getSecondOfPairFlag) fail(s"Read ${r2.getReadName} was not the second end of a pair")
        if (!r1.getReadName.equals(r2.getReadName)) fail(s"Read names did not match: '${r1.getReadName}' and '${r2.getReadName}'")

        // now do some work
        val bases1 = ExtractUmisFromBam.annotateRecord(record=r1, readStructure=rs1, molecularIndexTags=molecularIndexTags1, clippingAttribute=clippingAttribute)
        val bases2 = ExtractUmisFromBam.annotateRecord(record=r2, readStructure=rs2.get, molecularIndexTags=molecularIndexTags2, clippingAttribute=clippingAttribute)
        if (annotateReadNames) {
          // Developer Note: the delimiter must be an ascii character less than what is usually in the read names.  For
          // example "|" doesn't work.  I am not completely sure why.
          r1.setReadName(r1.getReadName + "+" + bases1 + bases2)
          r2.setReadName(r2.getReadName + "+" + bases1 + bases2)
        }
        assert(r1.getReadName.equals(r2.getReadName))

        // If we have duplicate tags, then concatenate them in the same order across the read pair.
        val tagAndValues = (molecularIndexTags1.map { tag => (tag, r1.getStringAttribute(tag)) }
          ++ molecularIndexTags2.map { tag => (tag, r2.getStringAttribute(tag)) }).toList
        tagAndValues.groupBy(_._1).foreach { case (tag, tuples) =>
          val attr = tuples.map(_._2).mkString(ExtractUmisFromBam.UmiDelimiter)
          r1.setAttribute(tag, attr)
          r2.setAttribute(tag, attr)
        }

        // If we have a single-tag, then also output values there
        singleTag.foreach { tag =>
          val value = tagAndValues.map { case (t,v) => v }.mkString(ExtractUmisFromBam.UmiDelimiter)
          r1.setAttribute(tag, value)
          r2.setAttribute(tag, value)
        }

        out.addAlignment(r1)
        out.addAlignment(r2)
        progress.record(r1, r2)

      }
      else {
        // verify everything is in order for single end reads
        if (!r1.getReadUnmappedFlag) fail(s"Read ${r1.getReadName} was not unmapped")

        // now do some work
        val bases1 = ExtractUmisFromBam.annotateRecord(record=r1, readStructure=rs1, molecularIndexTags=molecularIndexTags1, clippingAttribute=clippingAttribute)
        if (annotateReadNames) {
          // Developer Note: the delimiter must be an ascii character less than what is usually in the read names.  For
          // example "|" doesn't work.  I am not completely sure why.
          r1.setReadName(r1.getReadName + "+" + bases1)
        }

        // If we have duplicate tags, then concatenate them in the same order across the read
        val tagAndValues = molecularIndexTags1.map { tag => (tag, r1.getStringAttribute(tag)) }
        tagAndValues.groupBy(_._1).foreach { case (tag, tuples) =>
          val attr = tuples.map(_._2).mkString(ExtractUmisFromBam.UmiDelimiter)
          r1.setAttribute(tag, attr)
        }

        out.addAlignment(r1)
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
    * Extracts bases associated with molecular barcodes and adds them as SAM tags.  The read's bases will only contain
    * template bases as defined in the given read structure.
    */
  private[umi] def annotateRecord(record: SAMRecord,
                                  readStructure: ReadStructure,
                                  molecularIndexTags: Seq[String],
                                  clippingAttribute: Option[String] = None): String = {
    val bases = record.getReadString
    val qualities = record.getBaseQualityString
    // get the bases associated with each segment
    val readStructureBases = readStructure.structureReadWithQualities(bases, qualities, strict = false)
    // get the molecular barcode segments
    val molecularBarcodeBases = readStructureBases.collect { case SubRead(b: String, _, m: MolecularBarcode) => b }
    // set the barcode tags
    molecularIndexTags match {
      case Seq(tag) => record.setAttribute(tag, molecularBarcodeBases.mkString(UmiDelimiter))
      case _ =>
        if (molecularIndexTags.length < molecularBarcodeBases.length) throw new IllegalStateException("Found fewer molecular barcode SAM tags than molecular barcodes in the read structure.")
        else if (molecularIndexTags.length > molecularBarcodeBases.length) throw new IllegalStateException("Found fewer molecular barcodes in the read structure than molecular barcode SAM tags.")
        molecularIndexTags.zip(molecularBarcodeBases).foreach { case (tag, b) => record.setAttribute(tag, b) }
    }
    // keep only template bases and qualities in the output read
    val basesAndQualities = readStructureBases.flatMap { r =>
      r.segment match {
        case m: Template => Some((r.bases, r.quals.getOrElse(unreachable())))
        case _ => None
      }
    }
    // update any clipping information
    updateClippingInformation(record=record, clippingAttribute=clippingAttribute, readStructure=readStructure)
    record.setReadString(basesAndQualities.map(_._1).mkString)
    record.setBaseQualityString(basesAndQualities.map(_._2).mkString)
    // return the concatenation of the molecular barcodes in sequencing order
    molecularBarcodeBases.mkString
  }


  /**
    * Update the clipping information produced by Picard's MarkIlluminaAdapters to account for any non-template bases
    * that will be removed from the read.
    */
  private[umi] def updateClippingInformation(record: SAMRecord,
                                             clippingAttribute: Option[String],
                                             readStructure: ReadStructure): Unit = {
    clippingAttribute.map(tag => (tag, record.getIntegerAttribute(tag))) match {
      case None => Unit
      case Some((tag, null)) => Unit
      case Some((tag, clippingPosition)) =>
        val newClippingPosition = readStructure.takeWhile(_.offset < clippingPosition).collect { case t: Template =>
          if (t.offset + t.length < clippingPosition) t.length
          else clippingPosition - t.offset
        }.sum
        record.setAttribute(tag, newClippingPosition)
    }
  }
}
