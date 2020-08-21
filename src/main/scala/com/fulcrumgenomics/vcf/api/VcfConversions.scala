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

import java.util
import java.util.{List => JavaList}

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceMetadata
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api.VcfCount.Fixed
import htsjdk.variant.variantcontext.{GenotypeBuilder, VariantContext, VariantContextBuilder, Allele => JavaAllele}
import htsjdk.variant.vcf._

import scala.collection.JavaConverters.mapAsJavaMapConverter
import scala.collection.immutable.ListMap
import scala.collection.mutable

/**
  * Object that provides methods for converting from fgbio's scala VCF classes to HTSJDK's
  * Java VCF-related classes and vice-versa.  Only intended for internal use within the `api`
  * package and should not be exposed publicly.
  */
private[api] object VcfConversions {
  /** Class to allow wrapping an Array into an immutable IndexedSeq without copying. */
  private class ArrayIndexedSeq[T](private val array: Array[T]) extends scala.collection.immutable.IndexedSeq[T] {
    override final def apply(i: Int): T = this.array(i)
    override final def length: Int = this.array.length
    override final def nonEmpty: Boolean = this.array.length > 0
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
      val rec = c.getSAMSequenceRecord
      val length = if (rec.getSequenceLength == SequenceMetadata.UnknownSequenceLength) None else Some(rec.getSequenceLength)
      VcfContigHeader(rec.getSequenceIndex, rec.getSequenceName, length, Option(rec.getAssembly))
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
        // This is a horrible hack necessary because HTSJDK doesn't parse compound header lines unless
        // they are one of INFO/CONTIG/FILTER/ALT, and there's no way to get HTSJDK to report the version of the file!
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
      samples = in.getGenotypeSamples.toIndexedSeq
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
      val fields = new util.HashMap[String,String]()
      fields.put("ID", i.name)
      i.assembly.foreach(a => fields.put("assembly", a))
      i.length.foreach(l => fields.put("length", l.toString))
      out.addMetaDataLine(new VCFContigHeaderLine(fields, i.index))
    }

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
      case VCFHeaderLineCount.UNBOUNDED => VcfCount.Unknown
      case VCFHeaderLineCount.INTEGER   => VcfCount.Fixed(in.getCount)
    }
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
      case VcfCount.OnePerAltAllele    => Left(VCFHeaderLineCount.A)
      case VcfCount.OnePerGenotype     => Left(VCFHeaderLineCount.G)
      case VcfCount.OnePerAllele       => Left(VCFHeaderLineCount.R)
      case VcfCount.Unknown            => Left(VCFHeaderLineCount.UNBOUNDED)
      case VcfCount.Fixed(n)           => Right(n)
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
    * @return a [[Variant]] instance that is a copy of the [[VariantContext]] and does not rely on it
    *         post-return
    */
  def toScalaVariant(in: VariantContext, header: VcfHeader): Variant = try {
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

      val attrs = if (g.getExtendedAttributes.isEmpty && !g.hasAD && !g.hasDP && !g.hasGQ && !g.hasPL) Variant.EmptyGtAttrs else {
        val builder = Map.newBuilder[String, Any]
        if (g.hasAD) builder += ("AD" -> g.getAD.toIndexedSeq)
        if (g.hasDP) builder += ("DP" -> g.getDP)
        if (g.hasGQ) builder += ("GQ" -> g.getGQ)
        if (g.hasPL) builder += ("PL" -> g.getPL.toIndexedSeq)

        g.getExtendedAttributes.keySet().foreach { key =>
          val value = g.getExtendedAttribute(key)

          header.format.get(key) match {
            case Some(hd) => toTypedValue(value, hd.kind, hd.count).foreach(v => builder += (key -> v))
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
          case Some(hd) => toTypedValue(value, hd.kind, hd.count).foreach(v => builder += (key -> v))
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
    case ex: Throwable => throw new RuntimeException(s"Failed to convert variant: ${in}", ex)
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
    * @param value the value that came out of the INFO or genotype attributes from HTSJDK
    * @param kind the scala [[VcfFieldType]] declaring what the target type is
    * @param count the scala [[VcfCount]] describing how many values are expected
    * @return either a String/Float/Int/Char or an IndexedSeq of one of those types
    */
  private def toTypedValue(value: Any, kind: VcfFieldType, count: VcfCount): Option[Any] = ((value, kind, count): @unchecked) match {
    case (_, VcfFieldType.Flag, _       )              => Some(Variant.FlagValue)
    case (_, _,                 Fixed(0))              => Some(Variant.FlagValue)
    case (s: String, _,         Fixed(1))              => if (s == ".") None else Some(kind.parse(s))
    case (s: String, _,         _       )              => val xs = s.split(","); if (xs.forall(_ == ".")) None else Some(ArrayAttr(xs.map(kind.parse)))
    case (l: JavaList[String @unchecked], _, Fixed(1)) => if (l.get(0) == ".") None else Some(kind.parse(l.get(0)))
    case (l: JavaList[String @unchecked], _, _)        => if (l.forall(_ == ".")) None else Some(ArrayAttr(l.map(kind.parse)))
  }

  /**
    * Converts a Scala [[Variant]] back into a [[VariantContext]].
    *
    * @param in the [[Variant]] instance to convert
    * @param header the scala [[VcfHeader]] for the VCF being read or written
    * @return a VariantContext instance
    */
  def toJavaVariant(in: Variant, header: VcfHeader): VariantContext = {
    val alleles   = in.alleles.iterator.map { a => JavaAllele.create(a.toString, a eq in.alleles.ref) }.toJavaList
    val builder = new VariantContextBuilder(null, in.chrom, in.pos, in.end, alleles)
    in.id.foreach(i => builder.id(i))
    in.qual.foreach(q => builder.log10PError(q / -10 ))
    if (in.filters.isEmpty) builder.unfiltered() else builder.filters(in.filters.iterator.toJavaSet)
    builder.attributes(toJavaAttributeMap(in.attrs))

    val genotypes = header.samples.iterator.map { s =>
      val sgt = in.genotypes(s)
      val jgt = new GenotypeBuilder(s, sgt.callIndices.iterator.map(i => if (i == -1) JavaAllele.NO_CALL else alleles.get(i)).toJavaList)
      jgt.phased(sgt.phased)

      sgt.attrs.foreach {
        case ("GQ", value: Int)                    => jgt.GQ(value)
        case ("DP", value: Int)                    => jgt.DP(value)
        case ("AD", value: Seq[Int @unchecked])    => jgt.AD(value.toArray)
        case ("PL", value: Seq[Int @unchecked])    => jgt.PL(value.toArray)
        case ("FT", value: Seq[String @unchecked]) => value.foreach(f => jgt.filter(f))
        case (key,  value: Seq[Any])               => jgt.attribute(key, value.toArray)
        case (key,  value: Any)                    => jgt.attribute(key, value)
      }

      jgt.make()
    }.toJavaList

    builder.genotypes(genotypes).make()
  }

  /**
    * Converts a map of attributes from either the INFO field of a sample genotype into a java map that
    * HTSJDK can handle.  Largely this means translating to a [[java.util.Map]] and converting any scala
    * collection types back to arrays in the values.
    */
  private def toJavaAttributeMap(attrs: Map[String,Any]): java.util.LinkedHashMap[String,Any] = {
    val out = new util.LinkedHashMap[String,Any]()
    attrs.foreach {
      case (key, value: Seq[Any]) => out.put(key, value.toArray)
      case (key, value)           => out.put(key, value)
    }

    out
  }
}
