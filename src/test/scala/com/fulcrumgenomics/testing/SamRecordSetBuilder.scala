/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.testing

import java.nio.file.Files
import java.text.DecimalFormat
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.testing.SamRecordSetBuilder._
import dagr.commons.CommonsDef.{PathToBam, unreachable}
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.samtools._

import scala.collection.JavaConversions._

object SamRecordSetBuilder {
  sealed trait Strand { val isNegative: Boolean }
  object Plus  extends Strand { val isNegative = false }
  object Minus extends Strand { val isNegative = false }
}

/**
  * A scala wrapper to SAMRecordSetBuilder to make testing a little less unpleasant!
  */
class SamRecordSetBuilder(val readLength: Int=100,
                          val baseQuality: Int=30,
                          val sortOrder: SortOrder = SortOrder.unsorted,
                          val readGroupId: Option[String] = None
                         ) extends Iterable[SAMRecord] {
  private val builder = new SAMRecordSetBuilder(sortOrder != SortOrder.unsorted, sortOrder)
  readGroupId.foreach { id =>
    val readGroup = new SAMReadGroupRecord(id)
    readGroup.setSample("Sample")
    builder.setReadGroup(readGroup)
  }
  builder.setReadLength(readLength)
  builder.setUseNmFlag(false)

  private val counter = new AtomicLong(0)
  private val format  = new DecimalFormat("0000")
  private def nextName = format.format(counter.getAndIncrement())

  /** Adds a pair of reads to the file. */
  def addPair(name: String = nextName,
              contig: Int  = 0,
              start1: Int,
              start2: Int,
              record1Unmapped: Boolean = false,
              record2Unmapped: Boolean = false,
              cigar1: String = readLength + "M",
              cigar2: String = readLength + "M",
              mapq1: Int = 60,
              mapq2: Int = 60,
              strand1: Strand = Plus,
              strand2: Strand = Minus,
              baseQuality: Int = this.baseQuality,
              attrs: Map[String,AnyRef] = Map.empty
             ): Seq[SAMRecord] = {

    // Use the parent class to do the bulk of the work
    val recs = builder.addPair(name, contig, start1, start2, record1Unmapped, record2Unmapped, cigar1, cigar2,
      strand1.isNegative, strand2.isNegative, baseQuality).toList

    // Adjust the mapping qualities
    recs match {
      case r1 :: r2 :: Nil =>
        if (!r1.getReadUnmappedFlag) r1.setMappingQuality(mapq1)
        if (!r2.getReadUnmappedFlag) r2.setMappingQuality(mapq2)
        SamPairUtil.setMateInfo(r1, r2, true)
      case _ => unreachable("List should always contain a pair of reads!")
    }

    // Add any optional attributes
    for (attr <- attrs; rec <- recs) {
      rec.setAttribute(attr._1, attr._2)
    }

    recs
  }

  /** Adds a non-paired read. */
  def addFrag(name: String = nextName,
              contig: Int  = 0,
              start: Int,
              unmapped: Boolean = false,
              cigar: String = readLength + "M",
              mapq: Int = 60,
              strand: Strand = Plus,
              baseQuality: Int = this.baseQuality,
              attrs: Map[String,AnyRef] = Map.empty) : SAMRecord = {
    val rec = this.builder.addFrag(name, contig, start, strand.isNegative, unmapped, cigar, quals(readLength, baseQuality), baseQuality)

    // Adjust the mapping qualities
    if (!unmapped) rec.setMappingQuality(mapq)
    attrs.foreach(a => rec.setAttribute(a._1, a._2))

    rec
  }

  /** Generates a base quality string. */
  private def quals(length: Int, qual: Int) = new String(Array.fill(length)(SAMUtils.phredToFastq(qual)))

  /** Gets the SAMFileHeader that will be written. */
  def header = this.builder.getHeader

  /** Gets the SAMSequenceDictionary that will be written. */
  def dict   = this.builder.getHeader.getSequenceDictionary

  /** Returns an iterator over the records that have been built. */
  override def iterator = this.builder.iterator()

  /** Returns the number of records that have been built. */
  override def size = this.builder.size()

  /** Adds all the records from another builder to this one. */
  def += (other: SamRecordSetBuilder) : Unit = {
    other.iterator.foreach(this.builder.addRecord)
  }

  /** Writes the contents of the record set to the provided file path. */
  def write(path: PathToBam) : Unit = {
    val writer = new SAMFileWriterFactory()
      .setCreateIndex(sortOrder==SortOrder.coordinate)
      .makeWriter(builder.getHeader, true, path.toFile, null)
    iterator.foreach(writer.addAlignment)
    writer.close()
  }

  /** Writes the contents to a temporary file that will be deleted when the JVM exits. */
  def toTempFile(deleteOnExit: Boolean = true): PathToBam = {
    val path = Files.createTempFile("SamRecordSet.", ".bam")
    if (deleteOnExit) path.toFile.deleteOnExit()
    write(path)
    path
  }
}
