/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.testing

import java.nio.file.Files
import java.text.DecimalFormat
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.testing.SamBuilder._
import htsjdk.samtools._

import scala.collection.mutable.ArrayBuffer
import scala.util.Random

object SamBuilder {
  sealed trait Strand {
    val isNegative: Boolean
    override def toString(): String = this.getClass.getSimpleName.replaceFirst("[$].*$", "")
  }
  object Plus  extends Strand { val isNegative = false }
  object Minus extends Strand { val isNegative = true }
}

/** Class to create sets of [[SamRecord]]s for testing. */
class SamBuilder(val readLength: Int=100,
                 val baseQuality: Int=30,
                 val sort: Option[SamOrder] = None,
                 val readGroupId: Option[String] = None,
                 sd: Option[SAMSequenceDictionary] = None,
                 seed: Int = 42
                ) extends Iterable[SamRecord] {

  // Setup the header, sequence dictionary and read group
  val header = new SAMFileHeader()
  sort.getOrElse(SamOrder.Unsorted).applyTo(header)

  { // Build the default dictionary
    val dict: SAMSequenceDictionary = sd.getOrElse {
      val seqs = (Range.inclusive(1, 22) ++ Seq("X", "Y", "M")).map { chr => new SAMSequenceRecord("chr" + chr, 200e6.toInt) }
      new SAMSequenceDictionary(seqs.toIterator.toJavaList)
    }
    header.setSequenceDictionary(dict)
  }

  /** Shorter accessor for the sequence dictionary. */
  def dict: SAMSequenceDictionary = header.getSequenceDictionary

  val rg = new SAMReadGroupRecord(readGroupId.getOrElse("A"))
  rg.setSample("Sample")
  header.addReadGroup(rg)

  private val random  = new Random(seed)
  private val bases   = "ACGT".toCharArray
  private val counter = new AtomicLong(0)
  private val format  = new DecimalFormat("0000")

  /** The collection of records being accumulated. */
  private val records = ArrayBuffer[SamRecord]()

  /** The default string of qualities to use. */
  protected val defaultQuals: String = (33 + baseQuality).toChar.toString * readLength

  /** Generate a sequential name. */
  protected def nextName: String = format.format(counter.getAndIncrement())

  /** Generate a random sequence of bases of length readLength. */
  protected def randomBases: String = {
    val chs = new Array[Char](readLength)
    forloop (from=0, until=readLength) { i =>
      chs(i) = bases(random.nextInt(bases.length))
    }

    new String(chs)
  }


  /** Adds a pair of reads to the file. */
  def addPair(name: String = nextName,
              bases1 : String = randomBases,
              bases2: String = randomBases,
              quals1: String = defaultQuals,
              quals2: String = defaultQuals,
              contig: Int = 0,
              start1: Int = SAMRecord.NO_ALIGNMENT_START,
              start2: Int = SAMRecord.NO_ALIGNMENT_START,
              unmapped1: Boolean = false,
              unmapped2: Boolean = false,
              cigar1: String = readLength + "M",
              cigar2: String = readLength + "M",
              mapq1: Int = 60,
              mapq2: Int = 60,
              strand1: Strand = Plus,
              strand2: Strand = Minus,
              attrs: Map[String,Any] = Map.empty
             ): Seq[SamRecord] = {

    require(bases1.length == quals1.length, "bases1 and quals1 were different lengths.")
    require(bases2.length == quals2.length, "bases2 and quals2 were different lengths.")
    val cig1 = Cigar(cigar1)
    val cig2 = Cigar(cigar2)
    require(unmapped1 || bases1.length == cig1.lengthOnQuery, "bases1 doesn't agree with cigar on length.")
    require(unmapped2 || bases2.length == cig2.lengthOnQuery, "bases2 doesn't agree with cigar on length.")

    val r1            = SamRecord(header)
    r1.pf             = true
    r1.name           = name
    r1.bases          = bases1
    r1.quals          = quals1
    r1.refIndex       = if (start1 == SAMRecord.NO_ALIGNMENT_START) SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX else contig
    r1.start          = start1
    r1.positiveStrand = strand1 == Plus
    r1.paired         = true
    r1.firstOfPair    = true
    r1.unmapped       = unmapped1
    if (!unmapped1) {
      r1.cigar = cig1
      r1.mapq  = mapq1
    }
    r1("RG") = this.rg.getId
    attrs.foreach { case (key, value) => r1(key) = value }

    val r2            = SamRecord(header)
    r2.pf             = true
    r2.name           = name
    r2.bases          = bases2
    r2.quals          = quals2
    r2.refIndex       = if (start2 == SAMRecord.NO_ALIGNMENT_START) SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX else contig
    r2.start          = start2
    r2.positiveStrand = strand2 == Plus
    r2.paired         = true
    r2.secondOfPair   = true
    r2.unmapped       = unmapped2
    if (!unmapped2) {
      r2.cigar = cig2
      r2.mapq  = mapq2
    }
    r2("RG") = this.rg.getId
    attrs.foreach { case (key, value) => r2(key) = value }

    SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
    val recs = Seq(r1, r2)
    this.records ++= recs
    recs
  }

  /** Adds a non-paired read. Returns an Option for convenience of calling map/foreach etc. */
  def addFrag(name: String = nextName,
              bases: String = randomBases,
              quals: String = defaultQuals,
              contig: Int = 0,
              start: Int  = SAMRecord.NO_ALIGNMENT_START,
              unmapped: Boolean = false,
              cigar: String = readLength + "M",
              mapq: Int = 60,
              strand: Strand = Plus,
              attrs: Map[String,Any] = Map.empty) : Option[SamRecord] = {

    val cig = Cigar.fromSam(cigar)
    require(bases.length == quals.length, "bases and quals must be same length.")
    require(unmapped || bases.length == cig.lengthOnQuery, "bases and cigar disagree on length")

    val r1            = SamRecord(header)
    r1.pf             = true
    r1.name           = name
    r1.bases          = bases
    r1.quals          = quals
    r1.refIndex       = if (start == SAMRecord.NO_ALIGNMENT_START) SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX else contig
    r1.start          = start
    r1.positiveStrand = strand == Plus
    r1.unmapped       = unmapped
    if (!unmapped) {
      r1.cigar = cig
      r1.mapq  = mapq
    }
    r1("RG") = this.rg.getId
    attrs.foreach { case (key, value) => r1(key) = value }

    this.records += r1
    Some(r1)
  }

  /** Returns an iterator over the records that have been built. */
  override def iterator: Iterator[SamRecord] = this.sort match {
    case None => this.records.iterator
    case Some(so) if so == SamOrder.Unsorted || so == SamOrder.Unknown => this.records.iterator
    case Some(so) => this.records.toIndexedSeq.sortBy(so.sortkey).iterator
  }

  /** Returns the number of records that have been built. */
  override def size: Int = this.records.size

  /** Adds all the records from another builder to this one. */
  def ++= (that: Iterable[SamRecord]) : Unit = this.records ++= that

  /** Writes the contents of the record set to the provided file path. */
  def write(path: PathToBam) : Unit = {
    val writer = SamWriter(path, header, sort=sort)
    writer ++= this.records
    writer.close()
  }

  /** Writes the contents to a temporary file that will be deleted when the JVM exits. */
  def toTempFile(deleteOnExit: Boolean = true): PathToBam = {
    val path = Files.createTempFile("SamRecordSet.", ".bam")
    if (deleteOnExit) path.toFile.deleteOnExit()
    write(path)
    path
  }

  /** Creates a SamReader over the records stored in a temporary file. */
  def toSource: SamSource = {
    val bam = toTempFile()
    SamSource(bam)
  }
}
