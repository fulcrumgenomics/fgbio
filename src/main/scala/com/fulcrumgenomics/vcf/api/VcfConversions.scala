/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.vcf.api

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceMetadata
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api.VcfCount.Fixed
import htsjdk.variant.variantcontext.{GenotypeBuilder, VariantContext, VariantContextBuilder, Allele => JavaAllele}
import htsjdk.variant.vcf._

import java.util
import java.util.{List => JavaList}
import scala.collection.immutable.ListMap
import scala.jdk.CollectionConverters.MapHasAsJava

/** Thrown when a value of a typed INFO or FORMAT field is not a valid value of its declared type (e.g. `1.5` or
  * `3000000000` for an Integer field), so the VCF is invalid.  The message names the type, the value, the field and the
  * record (and, when read with a [[VcfSource]], the VCF). */
private[api] final class InvalidVcfValueException(message: String, cause: Throwable = null)
  extends IllegalArgumentException(message, cause)

/**
  * Object that provides methods for converting from fgbio's scala VCF classes to HTSJDK's
  * Java VCF-related classes and vice-versa.  Only intended for internal use within the `api`
  * package and should not be exposed publicly.
  */
private[api] object VcfConversions {
  /** The key of the INFO field giving the end position of a variant. */
  private val EndKey: String = "END"

  /** Class to allow wrapping an Array into an immutable IndexedSeq without copying. */
  private class ArrayIndexedSeq[T](private val array: Array[T]) extends scala.collection.immutable.IndexedSeq[T] {
    override final def apply(i: Int): T = this.array(i)
    override final def length: Int = this.array.length
    override final def isEmpty: Boolean = this.array.length == 0
  }

  /** Converts a String into Option[String]. Returns `None` if the input string is either
    * `null`, the empty string or [[Variant.Missing]].
    */
  private def opt(value: String): Option[String] = {
    if (value == null || value.isEmpty || value == Variant.Missing) None else Some(value)
  }

  /** Converts a Java VCF header into a scala VCF header. */
  def toScalaHeader(in: VCFHeader): VcfHeader = {
    val contigs = in.getContigLines.map { c =>
      VcfContigHeader.fromSequenceMetadata(SequenceMetadata(c.getSAMSequenceRecord))
    }.toIndexedSeq

    val infos = in.getInfoHeaderLines.toIndexedSeq.sortBy(_.getID).map { i =>
      VcfInfoHeader(i.getID, toScalaCount(i), toScalaKind(i.getType), i.getDescription, opt(i.getSource), opt(i.getVersion))
    }

    val formats = in.getFormatHeaderLines.toIndexedSeq.sortBy(_.getID).map { f =>
      VcfFormatHeader(f.getID, toScalaCount(f), toScalaKind(f.getType), f.getDescription)
    }

    val others = in.getOtherHeaderLines.map {
      case line: VCFSimpleHeaderLine =>
        val attrs = line.getGenericFields.entrySet().filter(_.getKey != "ID").map(entry => entry.getKey -> entry.getValue).toMap
        VcfGeneralHeader(line.getKey, line.getID, attrs)
      case line: VCFHeaderLine if line.getValue.startsWith("<") && line.getValue.endsWith(">") =>
        // This is a horrible hack necessary because HTSJDK only parses compound header lines of the types it knows.
        // VCF 4 parsing is used whatever the file's version, as the VCF 3 parser needs the list of expected tags.
        val attrs = VCFHeaderLineTranslator.parseLine(VCFHeaderVersion.VCF4_2, line.getValue, util.Collections.emptyList())
        val scalaAttrs = attrs.entrySet().filter(_.getKey != "ID").map(entry => entry.getKey -> entry.getValue).toMap
        VcfGeneralHeader(line.getKey, attrs.get("ID"), scalaAttrs)
      case line: VCFHeaderLine =>
        VcfGeneralHeader(line.getKey, line.getValue, Map.empty)
    }.toIndexedSeq

    VcfHeader(
      contigs = contigs,
      infos   = infos,
      formats = formats,
      filters = in.getFilterLines.toIndexedSeq.sortBy(_.getID).map(f => VcfFilterHeader(f.getID, f.getDescription)),
      others   = others,
      samples = in.getGenotypeSamples.toIndexedSeq,
      version = Option(in.getVCFHeaderVersion).map(toScalaVersion)
    )
  }


  /** Converts the scala VCF header back into a Java VCF header. */
  def toJavaHeader(in: VcfHeader): VCFHeader = {
    val out = new VCFHeader(new java.util.HashSet[VCFHeaderLine](), in.samples.iterator.toJavaList)

    in.infos.foreach { i =>
      val j = toJavaCount(i.count) match {
        case Left(countType) => new VCFInfoHeaderLine(i.id, countType, toJavaKind(i.kind), i.description, i.source.orNull, i.version.orNull)
        case Right(intCount) => new VCFInfoHeaderLine(i.id, intCount,  toJavaKind(i.kind), i.description, i.source.orNull, i.version.orNull)
      }
      out.addMetaDataLine(j)
    }

    in.formats.foreach { i =>
      val j = toJavaCount(i.count) match {
        case Left(countType) => new VCFFormatHeaderLine(i.id, countType, toJavaKind(i.kind), i.description)
        case Right(intCount) => new VCFFormatHeaderLine(i.id, intCount,  toJavaKind(i.kind), i.description)
      }
      out.addMetaDataLine(j)
    }

    in.filters.foreach { i =>  out.addMetaDataLine(new VCFFilterHeaderLine(i.id, i.description)) }

    in.others.foreach { i =>
      val j = if (i.data.isEmpty) new VCFHeaderLine(i.headerType, i.id) else {
        new VCFSimpleHeaderLine(i.headerType, (i.data ++ Map("ID" -> i.id)).asJava)
      }
      out.addMetaDataLine(j)
    }

    in.contigs.foreach { i =>
      // Ordered as HTSJDK writes contig lines built from a sequence dictionary
      val fields = new util.LinkedHashMap[String,String]()
      fields.put("ID", i.name)
      i.length.foreach(l => fields.put("length", l.toString))
      i.assembly.foreach(a => fields.put("assembly", a))
      i.md5.foreach(m => fields.put("md5", m))
      i.species.foreach(s => fields.put("species", s))
      i.url.foreach(u => fields.put("URL", u))
      out.addMetaDataLine(new VCFContigHeaderLine(fields, i.index))
    }

    // Set after adding the lines, since HTSJDK takes the version from any fileformat line added.  Without a declared
    // version, HTSJDK writes VCF 4.2 and refuses lines that need a newer one, so declare the version they need.
    val required = out.getMetaDataInInputOrder.iterator().map(_.minimumVersion).maxByOption(_.ordinal)
    in.version.map(toJavaVersion)
      .orElse(required.filter(VCFHeaderVersion.DEFAULT_VERSION.isOlderThan))
      .foreach(out.setVCFHeaderVersion)

    out
  }

  /**
    * Converts the count for a header line from Java to Scala.  Slighly complicated because HTSJDK
    * classes represent this as both an int `count` _and_ a enum, whereas in scala we always represent
    * it as a single object/instance of a trait.
    *
    * @param in the INFO or FORMAT header line from HTSJDK
    * @return a [[VcfCount]] instance
    */
  def toScalaCount(in: VCFCompoundHeaderLine): VcfCount = {
    in.getCountType match {
      case VCFHeaderLineCount.A         => VcfCount.OnePerAltAllele
      case VCFHeaderLineCount.G         => VcfCount.OnePerGenotype
      case VCFHeaderLineCount.R         => VcfCount.OnePerAllele
      case VCFHeaderLineCount.P         => VcfCount.OnePerGtAllele
      case VCFHeaderLineCount.LA        => VcfCount.OnePerLocalAltAllele
      case VCFHeaderLineCount.LR        => VcfCount.OnePerLocalAllele
      case VCFHeaderLineCount.LG        => VcfCount.OnePerLocalGenotype
      case VCFHeaderLineCount.M         => VcfCount.OnePerBaseModification
      case VCFHeaderLineCount.UNBOUNDED => VcfCount.Unknown
      case VCFHeaderLineCount.INTEGER   => VcfCount.Fixed(in.getCount)
    }
  }

  /** Converts an HTSJDK VCF version to the scala [[VcfVersion]]. */
  def toScalaVersion(in: VCFHeaderVersion): VcfVersion = in match {
    case VCFHeaderVersion.VCF3_2 => VcfVersion.V3_2
    case VCFHeaderVersion.VCF3_3 => VcfVersion.V3_3
    case VCFHeaderVersion.VCF4_0 => VcfVersion.V4_0
    case VCFHeaderVersion.VCF4_1 => VcfVersion.V4_1
    case VCFHeaderVersion.VCF4_2 => VcfVersion.V4_2
    case VCFHeaderVersion.VCF4_3 => VcfVersion.V4_3
    case VCFHeaderVersion.VCF4_4 => VcfVersion.V4_4
    case VCFHeaderVersion.VCF4_5 => VcfVersion.V4_5
  }

  /** Converts the scala [[VcfVersion]] back into the HTSJDK equivalent. */
  def toJavaVersion(in: VcfVersion): VCFHeaderVersion = in match {
    case VcfVersion.V3_2 => VCFHeaderVersion.VCF3_2
    case VcfVersion.V3_3 => VCFHeaderVersion.VCF3_3
    case VcfVersion.V4_0 => VCFHeaderVersion.VCF4_0
    case VcfVersion.V4_1 => VCFHeaderVersion.VCF4_1
    case VcfVersion.V4_2 => VCFHeaderVersion.VCF4_2
    case VcfVersion.V4_3 => VCFHeaderVersion.VCF4_3
    case VcfVersion.V4_4 => VCFHeaderVersion.VCF4_4
    case VcfVersion.V4_5 => VCFHeaderVersion.VCF4_5
  }

  /**
    * Converts an HTSJDK enum representing the type of value in an INFO or FORMAT field to a scala enum.
    */
  def toScalaKind(in: VCFHeaderLineType): VcfFieldType = in match {
    case VCFHeaderLineType.Character => VcfFieldType.Character
    case VCFHeaderLineType.Flag      => VcfFieldType.Flag
    case VCFHeaderLineType.Float     => VcfFieldType.Float
    case VCFHeaderLineType.Integer   => VcfFieldType.Integer
    case VCFHeaderLineType.String    => VcfFieldType.String
  }

  /** Converts back from the scala [[VcfCount]] back to the information needed to specify the count
    * in HTSJDK.
    *
    * @param in the [[VcfCount]] representing how many values an INFO or FORMAT field should have
    * @return either a [[VCFHeaderLineCount]] when possible or an [[Int]] if the count is [[Fixed]]
    */
  def toJavaCount(in: VcfCount): Either[VCFHeaderLineCount, Int] = {
    in match {
      case VcfCount.OnePerAltAllele        => Left(VCFHeaderLineCount.A)
      case VcfCount.OnePerGenotype         => Left(VCFHeaderLineCount.G)
      case VcfCount.OnePerAllele           => Left(VCFHeaderLineCount.R)
      case VcfCount.OnePerGtAllele         => Left(VCFHeaderLineCount.P)
      case VcfCount.OnePerLocalAltAllele   => Left(VCFHeaderLineCount.LA)
      case VcfCount.OnePerLocalAllele      => Left(VCFHeaderLineCount.LR)
      case VcfCount.OnePerLocalGenotype    => Left(VCFHeaderLineCount.LG)
      case VcfCount.OnePerBaseModification => Left(VCFHeaderLineCount.M)
      case VcfCount.Unknown                => Left(VCFHeaderLineCount.UNBOUNDED)
      case VcfCount.Fixed(n)               => Right(n)
    }
  }

  /**
    * Converts back from the scala [[VcfFieldType]] into the java/HTSJDK equivalent.
    */
  def toJavaKind(in: VcfFieldType ): VCFHeaderLineType = in match {
    case VcfFieldType.Character => VCFHeaderLineType.Character
    case VcfFieldType.Flag      => VCFHeaderLineType.Flag
    case VcfFieldType.Float     => VCFHeaderLineType.Float
    case VcfFieldType.Integer   => VCFHeaderLineType.Integer
    case VcfFieldType.String    => VCFHeaderLineType.String
  }

  /**
    * Converts a [[VariantContext]] and all nested classes into a [[Variant]] and set of [[Genotype]]s.
    *
    * @param in the [[VariantContext]] to be converted
    * @param header the scala [[VcfHeader]] which contains the definitions of all the INFO and FORMAT
    *               fields as well as the ordered list of sample names.
    * @param allowMissingFieldsInHeader if true, INFO and FORMAT fields not described in the header are kept as
    *                                   [[VcfFieldType.String]] fields of [[VcfCount.Unknown]] count (or as
    *                                   [[Variant.FlagValue]] for an INFO flag, and as an Int for the INFO END field),
    *                                   otherwise an exception is thrown
    * @return a [[Variant]] instance that is a copy of the [[VariantContext]] and does not rely on it
    *         post-return
    */
  def toScalaVariant(in: VariantContext, header: VcfHeader, allowMissingFieldsInHeader: Boolean = false): Variant = try {
    // Build up the allele set
    val scalaAlleles = in.getAlleles.iterator.map(a => Allele(a.getDisplayString)).toIndexedSeq
    val alleleMap    = in.getAlleles.iterator().zip(scalaAlleles).toMap
    val alleles      = AlleleSet(scalaAlleles.head, scalaAlleles.tail)

    // Build up the genotypes
    val gts = new Array[Genotype](in.getNSamples)
    forloop (from=0, until=in.getNSamples) { sampleIndex =>
      val g = in.getGenotype(sampleIndex)

      val calls = {
        val buffer = new Array[Allele](g.getPloidy)
        val javaAlleles = g.getAlleles
        forloop (from=0, until=buffer.length) { alleleIndex =>
          val a = javaAlleles.get(alleleIndex)
          buffer(alleleIndex) = if (a.isNoCall) NoCallAllele else alleleMap(a)
        }

        new ArrayIndexedSeq(buffer)
      }

      val attrs = if (g.getExtendedAttributes.isEmpty && !g.hasAD && !g.hasDP && !g.hasGQ && !g.hasPL && !g.isFiltered) Variant.EmptyGtAttrs else {
        val builder = Map.newBuilder[String, Any]
        if (g.hasAD) builder += ("AD" -> g.getAD.toIndexedSeq)
        if (g.hasDP) builder += ("DP" -> g.getDP)
        if (g.hasGQ) builder += ("GQ" -> g.getGQ)
        if (g.hasPL) builder += ("PL" -> g.getPL.toIndexedSeq)
        // HTSJDK holds the genotype filters (FT) apart from the other attributes, as one `;`-separated String
        if (g.isFiltered) builder += ("FT" -> g.getFilters)

        g.getExtendedAttributes.keySet().foreach { key =>
          val value = g.getExtendedAttribute(key)

          header.format.get(key) match {
            case Some(hd) =>
              val typed = try toTypedValue(value, hd.kind, hd.count) catch { case ex: InvalidVcfValueException =>
                throw new InvalidVcfValueException(s"${ex.getMessage} for FORMAT field $key of sample ${g.getSampleName} at ${in.getContig}:${in.getStart}")
              }
              typed.foreach(v => builder += (key -> v))
            case None if allowMissingFieldsInHeader => toUndescribedValue(value).foreach(v => builder += (key -> v))
            case None     => throw new IllegalStateException(s"Format field $key not described in header.")
          }
        }

        builder.result()
      }

      gts(sampleIndex) = Genotype(alleles, g.getSampleName, calls, g.isPhased, attrs)
    }

    // Build up the variant
    val inInfo = in.getAttributes
    val info = if (inInfo.isEmpty) Variant.EmptyInfo else {
      val builder = ListMap.newBuilder[String, Any]
      inInfo.entrySet().foreach { entry =>
        val key   = entry.getKey
        val value = entry.getValue
        header.info.get(key) match {
          case Some(hd) =>
            val typed = try toTypedValue(value, hd.kind, hd.count) catch { case ex: InvalidVcfValueException =>
              throw new InvalidVcfValueException(s"${ex.getMessage} for INFO field $key at ${in.getContig}:${in.getStart}")
            }
            typed.foreach(v => builder += (key -> v))
          // NB: Variant reads END as an Int, so read an END not described in the header as one
          case None if allowMissingFieldsInHeader && key == EndKey =>
            val typed = try toTypedValue(value, VcfFieldType.Integer, Fixed(1)) catch { case ex: InvalidVcfValueException =>
              throw new InvalidVcfValueException(s"${ex.getMessage} for INFO field $key at ${in.getContig}:${in.getStart}")
            }
            typed.foreach(v => builder += (key -> v))
          case None if allowMissingFieldsInHeader => toUndescribedValue(value).foreach(v => builder += (key -> v))
          case None     => throw new IllegalStateException(s"INFO field $key not described in header.")
        }
      }

      builder.result()
    }

    val filters = {
      if (in.filtersWereApplied() && in.isNotFiltered) Variant.PassingFilters
      else if (!in.filtersWereApplied()) Variant.EmptyFilters
      else in.getFilters.toSet
    }

    Variant(
      chrom     = in.getContig,
      pos       = in.getStart,
      id        = Option(if (in.getID == Variant.Missing) null else in.getID),
      alleles   = alleles,
      qual      = if (in.hasLog10PError) Some(in.getPhredScaledQual) else None,
      filters   = filters,
      attrs     = info,
      genotypes = new GenotypeMap(gts, header.sampleIndex)
    )
  }
  catch {
    case ex: InvalidVcfValueException => throw ex
    case ex: Throwable                => throw new RuntimeException(s"Failed to convert variant: ${in}", ex)
  }

  /**
    * Converts a value found in an INFO or FORMAT/genotype attribute to a well typed value. Because HTSJDK can
    * return any number of less-than-useful types, including strings, comma-separated strings and lists of strings
    * this has to do a bunch of matching to figure out how to convert both to the appropriate atomic type (e.g.
    * Int or Float) _and_ whether to return a singular value or a collection.
    *
    * Collection/array values are always returned as [[IndexedSeq]]s since those allow for indexed access and,
    * unlike arrays, are immutable.
    *
    * Note that in the case of Flag types and things that have a fixed count of 0, [[Variant.FlagValue]] is returned
    * since the value itself is unimportant.
    *
    * HTSJDK gives the values read from a BCF as typed values (e.g. Integer, Double or a list of them, with `null` for a
    * missing value) rather than as text, so these are first turned back into text.
    *
    * @param value the value that came out of the INFO or genotype attributes from HTSJDK
    * @param kind the scala [[VcfFieldType]] declaring what the target type is
    * @param count the scala [[VcfCount]] describing how many values are expected
    * @return either a String/Float/Int/Char or an IndexedSeq of one of those types
    */
  private def toTypedValue(value: Any, kind: VcfFieldType, count: VcfCount): Option[Any] = ((toText(value), kind, count): @unchecked) match {
    case (_, VcfFieldType.Flag, _       )                       => Some(Variant.FlagValue)
    case (_, _,                 Fixed(0))                       => Some(Variant.FlagValue)
    case (s: String, _,         Fixed(1))                       => if (s == ".") None else Some(parseValue(s, kind))
    case (s: String, _,         _       )                       => val xs = s.split(","); if (xs.forall(_ == ".")) None else Some(ArrayAttr(xs.map(parseValue(_, kind))))
    // A single String value that contains commas is split into a list by HTSJDK, so join it back up
    case (l: JavaList[String @unchecked], VcfFieldType.String, Fixed(1)) => val s = l.mkString(","); if (s == ".") None else Some(s)
    case (l: JavaList[String @unchecked], _, Fixed(1))          => if (l.get(0) == ".") None else Some(parseValue(l.get(0), kind))
    case (l: JavaList[String @unchecked], _, _)                 => if (l.forall(_ == ".")) None else Some(ArrayAttr(l.map(parseValue(_, kind))))
  }

  /** Parses a single value of the given type, throwing an [[InvalidVcfValueException]] that names the type and the value
    * if it is not a valid value of the type (e.g. `1.5` or `3000000000` for an Integer field). */
  private def parseValue(s: String, kind: VcfFieldType): Any = {
    try kind.parse(s) catch {
      case _: NumberFormatException => throw new InvalidVcfValueException(s"Invalid $kind value '$s'")
    }
  }

  /** Converts a value from HTSJDK to a String, or a list of Strings, as it would appear in a VCF.  Values read from a VCF
    * are already Strings (or lists of them); values read from a BCF are typed, with `null` for missing values. */
  private def toText(value: Any): Any = value match {
    case s: String                => s
    case null                     => Variant.Missing
    case l: JavaList[_]           => l.iterator().map(toText).map(_.toString).toJavaList
    case a: Array[_]              => a.iterator.map(toText).map(_.toString).toJavaList
    case d: java.lang.Double      => toText(d.floatValue()) // BCF stores 32-bit floats, which HTSJDK widens to Double
    case f: java.lang.Float       => f.toString
    case other                    => other.toString
  }

  /**
    * Converts the value of an INFO or FORMAT field that is not described in the header.  HTSJDK gives `true` for a
    * flag, which becomes [[Variant.FlagValue]]; any other value is treated as a [[VcfFieldType.String]] field with
    * [[VcfCount.Unknown]] count.
    */
  private def toUndescribedValue(value: Any): Option[Any] = value match {
    case _: java.lang.Boolean => Some(Variant.FlagValue)
    case _                    => toTypedValue(value, VcfFieldType.String, VcfCount.Unknown)
  }

  /**
    * Converts a Scala [[Variant]] back into a [[VariantContext]].
    *
    * @param in the [[Variant]] instance to convert
    * @param header the scala [[VcfHeader]] for the VCF being read or written
    * @param isBcf true if the variant will be written to a BCF, which needs numbers rather than text for Float values
    * @return a VariantContext instance
    */
  def toJavaVariant(in: Variant, header: VcfHeader, isBcf: Boolean = false): VariantContext = {
    val alleles   = in.alleles.iterator.map { a => JavaAllele.create(a.toString, a eq in.alleles.ref) }.toJavaList
    val builder = new VariantContextBuilder(null, in.chrom, in.pos.toLong, in.end.toLong, alleles)
    in.id.foreach(i => builder.id(i))
    in.qual.foreach(q => builder.log10PError(q / -10 ))
    if (in.filters.isEmpty) builder.unfiltered() else builder.filters(in.filters.iterator.toJavaSet)
    builder.attributes(toJavaInfoMap(in.attrs, header, isBcf))

    val genotypes = header.samples.iterator.map { s =>
      val sgt = in.genotypes(s)
      val jgt = new GenotypeBuilder(s, sgt.callIndices.iterator.map(i => if (i == -1) JavaAllele.NO_CALL else alleles.get(i)).toJavaList)
      jgt.phased(sgt.phased)

      sgt.attrs.foreach {
        case ("GQ", value: Int)                    => if (value != Variant.MissingInt) jgt.GQ(value)
        case ("DP", value: Int)                    => if (value != Variant.MissingInt) jgt.DP(value)
        case ("AD", value: Seq[Int @unchecked])    => jgt.AD(value.toArray)
        case ("PL", value: Seq[Int @unchecked])    => jgt.PL(value.toArray)
        case ("FT", value: String)                 => jgt.filter(value)
        case ("FT", value: Seq[String @unchecked]) => jgt.filters(value.iterator.toJavaList)
        case (key,  value: Any)                    => jgt.attribute(key, toJavaValue(value, isBcf))
      }

      jgt.make()
    }.toJavaList

    builder.genotypes(genotypes).make()
  }

  /** Converts the INFO attributes of a variant into a java map that HTSJDK can handle (see [[toJavaValue]]). */
  private def toJavaInfoMap(attrs: Map[String,Any], header: VcfHeader, isBcf: Boolean): java.util.LinkedHashMap[String,Any] = {
    val out = new util.LinkedHashMap[String,Any]()
    attrs.foreach {
      // An INFO field not described in the header that holds the flag value was read as a flag, so write it as one
      case (key, Variant.FlagValue) if !header.info.contains(key) => out.put(key, java.lang.Boolean.TRUE)
      case (key, value)                                           => out.put(key, toJavaValue(value, isBcf))
    }
    out
  }

  /**
    * Converts an INFO or genotype attribute value into a value that HTSJDK writes correctly:
    *   - a missing value ([[Variant.MissingInt]], [[Variant.MissingFloat]] or [[Variant.Missing]]), on its own or in a
    *     collection, becomes `null`, which HTSJDK writes as `.`
    *   - a Float is written as text with [[formatFloat]] in a VCF, or as a Double in a BCF (which HTSJDK requires)
    *   - a collection becomes an array for a VCF, or a [[java.util.List]] for a BCF
    */
  private def toJavaValue(value: Any, isBcf: Boolean): Any = {
    def element(x: Any): Any = x match {
      case _ if Variant.isMissingValue(x) => null
      case f: Float                       => if (isBcf) f.toDouble else formatFloat(f)
      case other                          => other
    }

    value match {
      case xs: Seq[Any] => if (isBcf) xs.iterator.map(element).toJavaList else xs.iterator.map(element).toArray[Any]
      case other        => element(other)
    }
  }

  /**
    * Formats a Float as the shortest decimal text that reads back as the same Float, without an exponent and without
    * trailing zeros (e.g. `60` rather than `60.0`, and `0.0001` rather than `1.0E-4`).  NaN and infinities are written
    * as `NaN`, `Inf` and `-Inf`, as the VCF specification describes.
    */
  private[api] def formatFloat(f: Float): String = {
    if (f.isNaN) "NaN"
    else if (f.isPosInfinity) "Inf"
    else if (f.isNegInfinity) "-Inf"
    else if (f == 0f) { if (java.lang.Float.compare(f, 0f) < 0) "-0" else "0" } // keep the sign of negative zero
    else BigDecimal(java.lang.Float.toString(f)).bigDecimal.stripTrailingZeros.toPlainString
  }
}
