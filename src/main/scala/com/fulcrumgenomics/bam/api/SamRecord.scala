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

package com.fulcrumgenomics.bam.api


import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.alignment.Cigar
import htsjdk.samtools
import htsjdk.samtools.SamPairUtil.PairOrientation
import htsjdk.samtools._
import htsjdk.samtools.util.CoordMath

// TODO long-term: methods for 5'end, 3' end, unclipped 5' end and unclipped 3' end
// TODO long-term: replacement for alignment blocks?

/**
  * Trait that extends [[SAMRecord]] so that it can override methods and call super.  Is separate from
  * [[SamRecord]] so that [[SamRecord]] does not export the [[SAMRecord]] API.
  *
  * Currently this class is only necessary because a) we define our own [[Cigar]] class and cache the
  * parsed value in [[SamRecord]] and b) we use HTSJDK utilities that may set the cigar using the
  * [[SAMRecord]] API methods `setCigar()` and `setCigarString()`.
  */
private[api] trait SamRecordIntermediate extends SAMRecord {
  final override def setCigar(cigar: samtools.Cigar): Unit = {
    val s = if (cigar == null) null else TextCigarCodec.encode(cigar)
    setCigarString(s)
  }

  override def setCigarString(cigar: String): Unit = { super.setCigarString(cigar);  cigarChanged() }
  protected[api] def setCigarStringNoNotify(cigar: String): Unit = super.setCigarString(cigar)
  protected[api] def cigarChanged(): Unit = {}
}

/** Class that is used to provide a nice API to transient attributes in the SamRecord. */
class TransientAttrs(private val rec: SamRecord) {
  def apply[A](key: Any): A = rec.asSam.getTransientAttribute(key).asInstanceOf[A]
  def update(key: Any, value: Any): Any = if (value == null) rec.asSam.removeTransientAttribute(key) else rec.asSam.setTransientAttribute(key, value)
  def get[A](key: Any): Option[A] = Option[A](apply[A](key))
  def getOrElse[A](key: Any, default: => A): A = rec.asSam.getTransientAttribute(key) match {
    case null => default
    case value => value.asInstanceOf[A]
  }
  def getOrElseUpdate[A](key: Any, default: => A): A = rec.asSam.getTransientAttribute(key) match {
    case null  =>
      val value: A = default
      val _ = update(key, value)
      value
    case value => value.asInstanceOf[A]
  }
}

/**
  * A trait that fgbio uses as a replacement for [[htsjdk.samtools.SAMRecord]].  The trait is self-typed as a
  * `SamRecordIntermediate` which is a subclass of SAMRecord.  It is done this wasy so that
  * a) we can access superclass methods via `SamRecordIntermediate` but that self-typing here
  * instead of extending hides the [[htsjdk.samtools.SAMRecord]] API from users of the class.  The result is
  * always a [[htsjdk.samtools.SAMRecord]] but isn't seen as such without casting.
  */
trait SamRecord {
  this: SamRecordIntermediate =>

  private var _cigar: Cigar = _

  ////////////////////////////////////////////////////////////////////////////
  // Wrapper methods.
  ////////////////////////////////////////////////////////////////////////////

  @inline final def name: String = getReadName
  @inline final def name_=(name: String): Unit = setReadName(name)

  @inline final def flags: Int = getFlags

  @inline final def paired: Boolean = getReadPairedFlag
  @inline final def paired_=(paired: Boolean):Unit = setReadPairedFlag(paired)

  @inline final def unpaired: Boolean = !paired
  @inline final def unpaired_=(unpaired: Boolean): Unit = this.paired = !unpaired

  @inline final def properlyPaired: Boolean = getProperPairFlag
  @inline final def properlyPaired_=(paired: Boolean):Unit = setProperPairFlag(paired)

  @inline final def mapped: Boolean = !getReadUnmappedFlag
  @inline final def mapped_=(mapped: Boolean):Unit = setReadUnmappedFlag(!mapped)
  @inline final def unmapped: Boolean = getReadUnmappedFlag
  @inline final def unmapped_=(unmapped: Boolean):Unit = setReadUnmappedFlag(unmapped)

  @inline final def mateMapped: Boolean = !getMateUnmappedFlag
  @inline final def mateMapped_=(mapped: Boolean):Unit = setMateUnmappedFlag(!mapped)
  @inline final def mateUnmapped: Boolean = getMateUnmappedFlag
  @inline final def mateUnmapped_=(unmapped: Boolean):Unit = setMateUnmappedFlag(unmapped)

  @inline final def positiveStrand: Boolean = !getReadNegativeStrandFlag
  @inline final def positiveStrand_=(positive: Boolean):Unit = setReadNegativeStrandFlag(!positive)
  @inline final def negativeStrand: Boolean = getReadNegativeStrandFlag
  @inline final def negativeStrand_=(negative: Boolean):Unit = setReadNegativeStrandFlag(negative)

  @inline final def matePositiveStrand: Boolean = !getMateNegativeStrandFlag
  @inline final def matePositiveStrand_=(positive: Boolean):Unit = setMateNegativeStrandFlag(!positive)
  @inline final def mateNegativeStrand: Boolean = getMateNegativeStrandFlag
  @inline final def mateNegativeStrand_=(negative: Boolean):Unit = setMateNegativeStrandFlag(negative)

  @inline final def firstOfPair: Boolean = getFirstOfPairFlag
  @inline final def firstOfPair_=(first: Boolean):Unit = setFirstOfPairFlag(first)

  @inline final def secondOfPair: Boolean = getSecondOfPairFlag
  @inline final def secondOfPair_=(second: Boolean):Unit = setSecondOfPairFlag(second)

  @inline final def secondary: Boolean = isSecondaryAlignment
  @inline final def secondary_=(secondary: Boolean):Unit = setSecondaryAlignment(secondary)
  @inline final def primary: Boolean = !this.secondary
  @inline final def primary_=(primary: Boolean):Unit = this.secondary = !primary

  @inline final def pf: Boolean = !getReadFailsVendorQualityCheckFlag
  @inline final def pf_=(pf: Boolean):Unit = setReadFailsVendorQualityCheckFlag(!pf)

  @inline final def duplicate: Boolean = getDuplicateReadFlag
  @inline final def duplicate_=(dupe: Boolean):Unit = setDuplicateReadFlag(dupe)

  @inline final def supplementary: Boolean = getSupplementaryAlignmentFlag
  @inline final def supplementary_=(supplementary: Boolean):Unit = setSupplementaryAlignmentFlag(supplementary)

  @inline final def refName: String = getReferenceName
  @inline final def refName_=(name: String):Unit = setReferenceName(name)

  @inline final def refIndex: Int = getReferenceIndex
  @inline final def refIndex_=(index: Int):Unit = setReferenceIndex(index)

  @inline final def start: Int = getAlignmentStart
  @inline final def start_=(s: Int):Unit = setAlignmentStart(s)

  @inline final def end: Int = if (unmapped) SAMRecord.NO_ALIGNMENT_START else start + cigar.lengthOnTarget - 1

  @inline final def unclippedStart: Int = if (unmapped) SAMRecord.NO_ALIGNMENT_START else getUnclippedStart
  @inline final def unclippedEnd: Int   = if (unmapped) SAMRecord.NO_ALIGNMENT_START else getUnclippedEnd

  @inline final def unSoftClippedStart: Int = if (unmapped) SAMRecord.NO_ALIGNMENT_START else start - cigar.leadingSoftClippedBases
  @inline final def unSoftClippedEnd: Int   = if (unmapped) SAMRecord.NO_ALIGNMENT_START else end + cigar.trailingSoftClippedBases

  @inline final def mapq: Int = getMappingQuality
  @inline final def mapq_=(q: Int):Unit = setMappingQuality(q)

  // The cigar handling code here is rather ugly because we can't _both_ override setCigar/setCigarString
  // and call super.setCigar/setCigarString in a self-typed trait :(  So we need to leave those methods
  // alone, but also detect when they've been called and reset our cigar representation!
  @inline final def cigar: Cigar = { if (_cigar == null) _cigar = Cigar.fromSam(getCigarString);  _cigar }
  @inline final def cigar_=(cig: String): Unit = { this._cigar = null; setCigarStringNoNotify(cig) }
  @inline final def cigar_=(cig: Cigar): Unit  = { this._cigar = cig; setCigarStringNoNotify(if (cig.isEmpty) null else cig.toString()) }
  @inline final override def cigarChanged(): Unit = { this._cigar = null } // null out the cigar; lazily reconstruct if asked for again

  @inline final def mateRefName: String = getMateReferenceName
  @inline final def mateRefName_=(name: String):Unit = setMateReferenceName(name)

  @inline final def mateRefIndex: Int = getMateReferenceIndex
  @inline final def mateRefIndex_=(index: Int):Unit = setMateReferenceIndex(index)

  @inline final def mateStart: Int = getMateAlignmentStart
  @inline final def mateStart_=(s: Int):Unit = setMateAlignmentStart(s)

  @inline final def mateCigar: Option[Cigar] = {
    require(paired && mateMapped, "Cannot get a mate cigar on read without a mapped mate.")
    get[String]("MC").map(Cigar.apply)
  }
  @inline final def mateEnd: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate end position on read without a mapped mate.")
    get[String]("MC").map(cig => mateStart + Cigar(cig).lengthOnTarget - 1)
  }
  @inline final def mateUnclippedStart: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate unclipped start position on read without a mapped mate.")
    get[String]("MC").map { cig =>
      val cigar = Cigar(cig)
      mateStart - cigar.leadingClippedBases
    }
  }
  @inline final def mateUnclippedEnd: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate unclipped end position on read without a mapped mate.")
    get[String]("MC").map { cig =>
      val cigar = Cigar(cig)
      mateStart + cigar.lengthOnTarget - 1 + cigar.trailingClippedBases
    }
  }
  @inline final def mateUnSoftClippedStart: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate un-soft-clipped start position on read without a mapped mate.")
    get[String]("MC").map { cig =>
      val cigar = Cigar(cig)
      mateStart - cigar.leadingSoftClippedBases
    }
  }
  @inline final def mateUnSoftClippedEnd: Option[Int] = {
    require(paired && mateMapped, "Cannot get mate un-soft-clipped end position on read without a mapped mate.")
    get[String]("MC").map { cig =>
      val cigar = Cigar(cig)
      mateStart + cigar.lengthOnTarget - 1 + cigar.trailingSoftClippedBases
    }
  }
  @inline final def matesOverlap: Option[Boolean] = {
    require(mapped && paired && mateMapped, "Cannot determine if mates overlap without paired mates that are both mapped.")
    if (refIndex != mateRefIndex) Some(false)
    else if (mateStart > end) Some(false)
    else if (mateStart >= start && mateStart <= end) Some(true)
    else mateEnd.map(mateEnd => CoordMath.overlaps(start, end, mateStart, mateEnd))
  }

  @inline final def insertSize: Int = getInferredInsertSize
  @inline final def insertSize_=(s: Int):Unit = setInferredInsertSize(s)

  @inline final def bases: Array[Byte] = getReadBases
  @inline final def basesString: String = getReadString
  @inline final def bases_=(bs: Array[Byte]): Unit = setReadBases(bs)
  @inline final def bases_=(bs: String): Unit = setReadString(bs)

  @inline final def quals: Array[Byte] = getBaseQualities
  @inline final def qualsString: String = getBaseQualityString
  @inline final def quals_=(qs: Array[Byte]): Unit = setBaseQualities(qs)
  @inline final def quals_=(qs: String): Unit = setBaseQualityString(qs)

  @inline final def length: Int = getReadLength

  @inline final def header: SAMFileHeader = getHeader
  @inline final def header_=(header: SAMFileHeader) = setHeader(header)
  @inline final def readGroup: SAMReadGroupRecord = getReadGroup

  // Use apply/update for tag attributes
  @inline final def apply[A](name: String): A                    = getAttribute(name).asInstanceOf[A]
  @inline final def get[A](name: String): Option[A]              = Option[A](apply[A](name))
  @inline final def update(name: String, value: Any): Unit       = setAttribute(name, value)
  @inline final def attributes: Map[String,Any]                  = getAttributes.map(x => x.tag -> x.value).toMap
  @inline final def getOrElse[A](name: String, default: => A): A = get(name).getOrElse(default)
  @inline final def contains(name: String): Boolean              = hasAttribute(name)
  @inline final def remove(name: String): Unit                   = setAttribute(name, null)

  // transient attributes
  @inline final val transientAttrs: TransientAttrs = new TransientAttrs(this)

  // TODO long-term: replace these two methods with methods on [[Cigar]] to save creating alignment blocks in memory
  @inline final def refPosAtReadPos(pos: Int) = getReferencePositionAtReadPosition(pos)
  @inline final def readPosAtRefPos(pos: Int, returnLastBaseIfDeleted: Boolean) = getReadPositionAtReferencePosition(pos, returnLastBaseIfDeleted)

  /** Reconstructs the reference bases spanned by this record from its read bases, CIGAR, and MD tag.
    *
    * The returned bases span the alignment from its start to its end, inclusive: soft-clipped and inserted read bases
    * are omitted, while deleted reference bases are included, as are skipped (`N` CIGAR operator) reference bases as `N`s
    * since the MD tag does not store them.  Any letter is accepted for a mismatched or deleted reference base (e.g.
    * IUPAC codes).  Adjacent deletions in the CIGAR may be given as one run of deleted bases (e.g. `1^AC1`) or as
    * separate runs (e.g. `1^A0^C1`, as written by htsjdk and samtools).
    *
    * The record must not use `=` for a read base that the MD tag reports as matching the reference, since the
    * reference base cannot then be determined.
    *
    * @throws java.lang.IllegalArgumentException if the MD tag is missing or malformed, does not match the CIGAR, or the
    *                                            reference base for an `=` read base is needed
    */
  final def refBasesFromMd: Array[Byte] = {
    val md        = get[String]("MD").getOrElse(
      throw new IllegalArgumentException(s"Cannot get the reference bases without the MD tag for record with name: $name")
    )
    val readBases = this.bases
    val cigar     = this.cigar
    require(readBases.length == cigar.lengthOnQuery,
      s"Cannot get the reference bases from the MD tag without the read bases for record with name: $name")
    def mdFail(reason: String): Nothing =
      throw new IllegalArgumentException(s"MD tag '$md' does not match CIGAR '$cigar' ($reason) for record with name: $name")

    val refBases     = Array.newBuilder[Byte]
    var mdIndex      = 0     // the index of the next character to parse in the MD tag
    var matchesLeft  = 0     // the number of matching bases left from the last number parsed in the MD tag
    var inDeletion   = false // true when the MD tag is within a run of deleted bases (after a '^')
    var expectNumber = false // true after a mismatched or deleted base, which must be followed by a number
    var readOffset   = 0

    // True if the next character in the MD tag satisfies the predicate
    def nextIs(p: Char => Boolean): Boolean = mdIndex < md.length && p(md.charAt(mdIndex))

    // Parses any numbers in the MD tag until there are matching bases left or the next character is not a digit
    def parseMatches(): Unit = {
      while (matchesLeft == 0 && nextIs(Character.isDigit)) {
        var number = 0
        while (nextIs(Character.isDigit)) {
          number = number * 10 + (md.charAt(mdIndex) - '0')
          mdIndex += 1
        }
        matchesLeft  = number
        expectNumber = false
      }
    }

    cigar.foreach { elem =>
      if (elem.operator != CigarOperator.D) inDeletion = false
      elem.operator match {
        case CigarOperator.M | CigarOperator.EQ | CigarOperator.X =>
          forloop(from=0, until=elem.length) { i =>
            parseMatches()
            val readBase = readBases(readOffset + i)
            if (matchesLeft > 0) {
              require(readBase != '=', s"Cannot get the reference base from the MD tag for an '=' read base for record with name: $name")
              refBases += readBase
              matchesLeft -= 1
            }
            else if (nextIs(Character.isLetter)) {
              if (expectNumber) mdFail("expected a number between mismatched or deleted bases")
              refBases += md.charAt(mdIndex).toByte
              mdIndex     += 1
              expectNumber = true
            }
            else mdFail("expected a match or mismatch")
          }
        case CigarOperator.D =>
          // A deletion directly after another deletion may continue its run of deleted bases, or start a new run after a
          // zero (e.g. `^A0^C`).  Otherwise, a new run of deleted bases (after a '^') is required.
          if (!inDeletion || !nextIs(Character.isLetter)) {
            parseMatches()
            if (matchesLeft > 0 || !nextIs(_ == '^')) mdFail("expected a deletion")
            if (expectNumber && !inDeletion) mdFail("expected a number between mismatched or deleted bases")
            mdIndex   += 1
            inDeletion = true
          }
          forloop(from=0, until=elem.length) { _ =>
            if (!nextIs(Character.isLetter)) mdFail("expected a deleted base")
            refBases += md.charAt(mdIndex).toByte
            mdIndex += 1
          }
          expectNumber = true
        case CigarOperator.N =>
          // The MD tag does not store skipped bases
          forloop(from=0, until=elem.length) { _ => refBases += 'N'.toByte }
        case _ => () // I, S, H, and P do not consume reference bases
      }
      readOffset += elem.lengthOnQuery
    }

    parseMatches()
    if (matchesLeft > 0 || mdIndex < md.length) mdFail("the MD tag is longer than the alignment")
    refBases.result()
  }

  ////////////////////////////////////////////////////////////////////////////
  // Non-wrapper methods.
  ////////////////////////////////////////////////////////////////////////////

  /** Returns a string that is useful to identify a SamRecord, mostly for testing and error messages. */
  def id: String = {
    val builder = new StringBuilder
    builder.append(name)
    builder.append(if (paired && secondOfPair) "/2" else "/1")
    if (secondary) builder.append(":sec")
    if (supplementary) builder.append(":sup")
    builder.toString()
  }

  /** Returns this record as a SAMRecord. */
  def asSam: SAMRecord = this.asInstanceOf[SAMRecord]

  /** Gets the PairOrientation of the record. */
  def pairOrientation: PairOrientation = SamPairUtil.getPairOrientation(this.asSam)

  /** Returns true if the read is mapped in an FR pair, false otherwise. */
  def isFrPair: Boolean = {
    paired &&
      mapped &&
      mateMapped &&
      refIndex == mateRefIndex &&
      SamPairUtil.getPairOrientation(this) == PairOrientation.FR
  }

  /** Clone method that does a "reasonably deep" clone. The bases and quals are cloned as is the attributes map,
    * though not the values in the attributes map. */
  override def clone(): SamRecord = {
    val r = super.clone().asInstanceOf[SamRecord]
    r.bases = this.bases.clone()
    r.quals = this.quals.clone()
    r
  }
}

object SamRecord {

  /** The mapping quality when not available. */
  val UnknownMappingQuality: Int = SAMRecord.UNKNOWN_MAPPING_QUALITY

  /** The mapping quality for an unmapped or multi-mapped read */
  val ZeroMappingQuality: Int = SAMRecord.NO_MAPPING_QUALITY

  /** The reference name for an unmapped read. */
  val UnmappedReferenceName: String = SAMRecord.NO_ALIGNMENT_REFERENCE_NAME

  /** The reference index for an unmapped read. */
  val UnmappedReferenceIndex: Int = SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX

  /** The cigar string for an unmapped read. */
  val UnmappedCigarString: String = SAMRecord.NO_ALIGNMENT_CIGAR

  /** The reference position for an unmapped read. */
  val UnmappedStart: Int = SAMRecord.NO_ALIGNMENT_START

  /** The value for the base string when no bases are present. */
  val MissingBases: String = SAMRecord.NULL_SEQUENCE_STRING

  /** The value for the base quality string when no base qualities are present. */
  val MissingQuals: String = SAMRecord.NULL_QUALS_STRING

  /** The maximum insert size that can be stored in a [[SamRecord]]. */
  val MaximumInsertSize: Int = SAMRecord.MAX_INSERT_SIZE

  /** SAMRecord with mixin to add behaviour. */
  private final class EnhancedSamRecord(header: SAMFileHeader) extends SAMRecord(header) with SamRecordIntermediate with SamRecord

  /** BAMRecord with mixin to add behaviour. */
  private final class EnhancedBamRecord(header: SAMFileHeader,
                                        refIndex: Int,
                                        alignmentStart: Int,
                                        readNameLength: Short,
                                        mapq: Short,
                                        indexingBin: Int,
                                        cigarLen: Int,
                                        flags: Int,
                                        readLen: Int,
                                        mateRefIndex: Int,
                                        mateAlignmentStart: Int,
                                        insertSize: Int,
                                        variableLengthBlock: Array[Byte]) extends
    BAMRecord(header, refIndex, alignmentStart, readNameLength, mapq, indexingBin, cigarLen,
      flags, readLen, mateRefIndex, mateAlignmentStart, insertSize, variableLengthBlock) with SamRecordIntermediate with SamRecord

  /** Factory method for new [[SamRecord]] instances. */
  def apply(header: SAMFileHeader): SamRecord = new EnhancedSamRecord(header)

  /**
    * Singleton implementation of SAMRecordFactory to return instances of [[htsjdk.samtools.SAMRecord]] and
    * [[htsjdk.samtools.BAMRecord]] that have been enhanced by mixing in [[SamRecord]].
    */
  object Factory extends SAMRecordFactory {
    /** Factory implementation method to create SAMRecords. */
    override def createSAMRecord(header: SAMFileHeader): SAMRecord = new EnhancedSamRecord(header)

    /** Factory implementation method to create BAMRecords. */
    override def createBAMRecord(header: SAMFileHeader,
                                 refIndex: Int,
                                 alignmentStart: Int,
                                 readNameLength: Short,
                                 mapq: Short,
                                 indexingBin: Int,
                                 cigarLen: Int,
                                 flags: Int,
                                 readLen: Int,
                                 mateRefIndex: Int,
                                 mateAlignmentStart: Int,
                                 insertSize: Int,
                                 variableLengthBlock: Array[Byte]): BAMRecord =
      new EnhancedBamRecord(header, refIndex, alignmentStart, readNameLength, mapq, indexingBin, cigarLen,
        flags, readLen, mateRefIndex, mateAlignmentStart, insertSize, variableLengthBlock)
  }
}
