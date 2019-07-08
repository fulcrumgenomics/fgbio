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
import com.fulcrumgenomics.vcf.api.VcfCount.Fixed
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.variant.variantcontext.{GenotypeBuilder, VariantContext, VariantContextBuilder, Allele => JavaAllele}
import htsjdk.variant.vcf._

import scala.collection.JavaConverters.mapAsJavaMap
import scala.collection.immutable.ListMap
import scala.collection.mutable.ArrayBuffer

/**
  * Object that provides methods for converting from fgbio's scala VCF classes to HTSJDK's
  * Java VCF-related classes and vice-versa.
  */
private[api] object VcfConversions {
  /** Value used in VCF for values that are missing. */
  val Missing: String = "."

  /** Converts a String into Option[String] accounting for various empty/missing values. */
  private def opt(value: String): Option[String] = {
    if (value == null || value.isEmpty || value == Missing) None else Some(value)
  }

  /** Converts a Java VCF header into a scala VCF header. */
  def toScalaHeader(in: VCFHeader): VcfHeader = {
    val contigs = in.getContigLines.map { c =>
      val rec = c.getSAMSequenceRecord
      val length = if (rec.getSequenceLength == SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH) None else Some(rec.getSequenceLength)
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
        val attrs = line.getGenericFields.entrySet().map { entry => entry.getKey -> entry.getValue }.toMap
        VcfGeneralHeader(line.getKey, line.getID, attrs)
      case line: VCFHeaderLine =>
        VcfGeneralHeader(line.getKey, line.getValue, Map.empty)
    }.toIndexedSeq

    VcfHeader(
      contigs = contigs,
      infos   = infos,
      formats = formats,
      filters = in.getFilterLines.toIndexedSeq.sortBy(_.getID).map(f => VcfFilterHeader(f.getID, f.getDescription)),
      other   = others,
      samples = in.getSampleNamesInOrder.toIndexedSeq
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

    in.other.foreach { i =>
      val j = if (i.data.isEmpty) new VCFHeaderLine(i.headerType, i.id) else {
        new VCFSimpleHeaderLine(i.headerType, mapAsJavaMap(i.data ++ Map("ID" -> i.id)))
      }
      out.addMetaDataLine(j)
    }

    in.contigs.foreach { i =>
      val fields = new util.HashMap[String,String]()
      fields.put("ID", i.name)
      i.assembly.foreach(a => fields.put("assembly", a))
      out.addMetaDataLine(new VCFContigHeaderLine(fields, i.index))
    }

    out
  }

  def toScalaCount(in: VCFCompoundHeaderLine): VcfCount = {
    in.getCountType match {
      case VCFHeaderLineCount.A         => VcfCount.OnePerAltAllele
      case VCFHeaderLineCount.G         => VcfCount.OnePerGenotype
      case VCFHeaderLineCount.R         => VcfCount.OnePerAllele
      case VCFHeaderLineCount.UNBOUNDED => VcfCount.Unknown
      case VCFHeaderLineCount.INTEGER   => VcfCount.Fixed(in.getCount)
    }
  }

  def toScalaKind(in: VCFHeaderLineType): VcfFieldType = in match {
    case VCFHeaderLineType.Character => VcfFieldType.Character
    case VCFHeaderLineType.Flag      => VcfFieldType.Flag
    case VCFHeaderLineType.Float     => VcfFieldType.Float
    case VCFHeaderLineType.Integer   => VcfFieldType.Integer
    case VCFHeaderLineType.String    => VcfFieldType.String
  }

  def toJavaCount(in: VcfCount): Either[VCFHeaderLineCount, Int] = {
    in match {
      case VcfCount.OnePerAltAllele    => Left(VCFHeaderLineCount.A)
      case VcfCount.OnePerGenotype     => Left(VCFHeaderLineCount.G)
      case VcfCount.OnePerAllele       => Left(VCFHeaderLineCount.R)
      case VcfCount.Unknown            => Left(VCFHeaderLineCount.UNBOUNDED)
      case VcfCount.Fixed(n)           => Right(n)
    }
  }

  def toJavaKind(in: VcfFieldType ): VCFHeaderLineType= in match {
    case VcfFieldType.Character => VCFHeaderLineType.Character
    case VcfFieldType.Flag      => VCFHeaderLineType.Flag
    case VcfFieldType.Float     => VCFHeaderLineType.Float
    case VcfFieldType.Integer   => VCFHeaderLineType.Integer
    case VcfFieldType.String    => VCFHeaderLineType.String
  }

  def toScalaVariant(in: VariantContext, header: VcfHeader): Variant = {
    // Build up the allele set
    val ref  = Allele(in.getReference.getDisplayString)
    val alts = in.getAlternateAlleles.map(a => Allele(a.getDisplayString)).toIndexedSeq
    val alleles = AlleleSet(ref, alts)
    val alleleMap = alleles.map(a => a.toString -> a).toMap

    // Build up the genotypes
    val gts = in.getGenotypes.map { g =>
      val calls = g.getAlleles.map(a => alleleMap(a.getDisplayString)).toIndexedSeq
      val attrs = {
        val buffer = new ArrayBuffer[(String,Any)](g.getExtendedAttributes.size() + 4)
        if (g.hasAD) buffer.append("AD" -> g.getAD.toIndexedSeq)
        if (g.hasDP) buffer.append("DP" -> g.getDP)
        if (g.hasGQ) buffer.append("GQ" -> g.getGQ)
        if (g.hasPL) buffer.append("PL" -> g.getPL.toIndexedSeq)

        g.getExtendedAttributes.keySet().foreach { key =>
          val value = g.getExtendedAttribute(key)

          header.format.get(key) match {
            case Some(hd) => buffer.append(key -> toTypedValue(value, hd.kind, hd.count))
            case None     => throw new IllegalStateException(s"Format field $key not described in header.")
          }
        }

        buffer.toMap
      }

      Genotype(alleles, g.getSampleName, calls, g.isPhased, attrs)
    }.toIndexedSeq

    // Build up the variant
    val info = in.getAttributes.entrySet().map { entry =>
      val key   = entry.getKey
      val value = entry.getValue
      header.info.get(key) match {
        case Some(hd) => key -> toTypedValue(value, hd.kind, hd.count)
        case None     => throw new IllegalStateException(s"INFO field $key not described in header.")
      }
    }.toSeq

    Variant(
      chrom     = in.getContig,
      pos       = in.getStart,
      id        = Option(if (in.getID == Missing) null else in.getID),
      alleles   = alleles,
      qual      = if (in.hasLog10PError) Some(in.getPhredScaledQual) else None,
      filter    = in.getFilters.toIndexedSeq,
      info      = ListMap(info:_*),
      genotypes = gts.iterator.map(g => g.sample -> g).toMap
    )
  }

  def toTypedValue(value: Any, kind: VcfFieldType, count: VcfCount): Any = (value, kind, count) match {
    case (_, VcfFieldType.Flag, _       )   => None
    case (_, _,                 Fixed(0))   => None
    case (s: String, _,         Fixed(1))   => kind.parse(s)
    case (s: String, _,         _       )   => s.split(',').map(kind.parse).toIndexedSeq
    case (l: JavaList[String], _, Fixed(1)) => kind.parse(l.get(0))
    case (l: JavaList[String], _, _)        => l.map(kind.parse).toIndexedSeq
  }

  def toJavaVariant(in: Variant, header: VcfHeader): VariantContext = {
    val alleles   = in.alleles.iterator.map { a => JavaAllele.create(a.toString, a eq in.alleles.ref) }.toJavaList
    val builder = new VariantContextBuilder(null, in.chrom, in.pos, in.end, alleles)
    in.id.foreach(i => builder.id(i))
    in.qual.foreach(q => builder.log10PError(q / -10 ))
    if (in.filter.isEmpty) builder.unfiltered() else builder.filters(in.filter.iterator.toJavaSet)
    builder.attributes(toJavaAttributeMap(in.info))

    val genotypes = header.samples.iterator.map { s =>
      val sgt = in.genotypes(s)
      val jgt = new GenotypeBuilder(s, sgt.callIndices.iterator.map(i => alleles.get(i)).toJavaList)
      jgt.phased(sgt.phased)

      sgt.attributes.foreach {
        case ("GQ", value: Int)         => jgt.GQ(value)
        case ("DP", value: Int)         => jgt.DP(value)
        case ("AD", value: Seq[Int])    => jgt.AD(value.toArray)
        case ("PL", value: Seq[Int])    => jgt.PL(value.toArray)
        case ("FT", value: Seq[String]) => value.foreach(f => jgt.filter(f))
        case (key, value: Seq[Any])     => jgt.attribute(key, value.toArray)
        case (key, value: Any)          => jgt.attribute(key, value)
      }

      jgt.make()
    }.toJavaList

    builder.genotypes(genotypes).make()
  }

  def toJavaAttributeMap(attrs: Map[String,Any]): java.util.LinkedHashMap[String,Any] = {
    val out = new util.LinkedHashMap[String,Any]()
    attrs.foreach {
      case (key, value: Seq[Any]) => out.put(key, value.toArray)
      case (key, value)           => out.put(key, value)
    }

    out
  }
}
