/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.vcf.parsing.header

import com.fulcrumgenomics.vcf.api.VcfCount._
import com.fulcrumgenomics.vcf.api.{VcfAlternateAlleleHeader, VcfContigHeader, VcfCountAndKindHeader, VcfFilterHeader, VcfFormatHeader, VcfGeneralHeader, VcfHeaderEntry, VcfInfoHeader, VcfSymbolicAllele, VcfCount => Count, VcfFieldType => Kind}
import com.fulcrumgenomics.vcf.parsing.util.Lookup
import com.fulcrumgenomics.vcf.parsing.util.ParseResult.fail
import com.fulcrumgenomics.vcf.parsing.util.Parsing.csvKeyValues
import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex.SI
import fastparse.NoWhitespace._
import fastparse.{CharIn, CharPred, CharsWhile, CharsWhileIn, Fail, P, Start, StringIn, _}

object VcfHeaderParser {
  import MetaHeaderKey.{metaHeaderKeyToString, _}

  /** Parses a VCF Header entry/line. */
  def header[_: P]: P[VcfHeaderEntry] = Start ~ "##" ~ CharsWhile(_ != '=').!.flatMap {
    case "fileformat" => VcfFileFormatHeaderParser.parse
    case "contig"     => VcfContigHeaderParser.parse
    case "FILTER"     => VcfFilterHeaderParser.parse
    case "INFO"       => VcfInfoHeaderParser.parse
    case "FORMAT"     => VcfFormatHeaderParser.parse
    case "ALT"        => VcfAlternateAlleleHeaderParser.parse
    case headerType   => VcfGeneralHeaderParser.parse(headerType=headerType)
  }

  /** Parses the header line after the VCF Header returning zero or more sample names. */
  def samples[_: P]: P[Seq[String]] = {
    Start ~ "#CHROM" ~/ "\t" ~/ "POS" ~/ "\t" ~/ "ID" ~/ "\t" ~/ "REF" ~/ "\t" ~/ "ALT" ~/
      "\t" ~/ "QUAL" ~/ "\t" ~/ "FILTER" ~/ "\t" ~/ "INFO" ~/ "\t" ~/ "FORMAT" ~/ ("\t" ~/
      CharPred(_ != '\t').rep(1).!).rep ~ End
  }

  /** Parses the file format VCF header entry/line. */
  object VcfFileFormatHeaderParser {
    /** The valid versions supported by the parser. */
    def versions[_: P]: P[String] = StringIn("VCFv4.3").!
    /** Parses the remaining input after the leading header key. */
    def parse[_: P]: P[VcfGeneralHeader] = "=" ~/ versions.map { version =>
      VcfGeneralHeader(headerType="fileformat", id=version)
    } ~ End
  }

  /** Parses the [[VcfContigHeader]] entry/line. */
  object VcfContigHeaderParser {
    /** Parses the leading character in a contig name. */
    def leadingChar[_: P]: P[String] = (CharIn("0-9A-Za-z!#$%&+./:;?@^_|~") | CharPred(_ == '-')).!
    /** Parses zero or more trailing characters in a contig name. */
    def trailingString[_: P]: P[String] = (CharsWhileIn("0-9A-Za-z!#$%&*+./:;=?@^_|~", 0) | CharPred(_ == '-').rep(0)).!
    /** Parses a contig name. */
    def name[_: P]: P[String] = (leadingChar ~ trailingString).map {
      case (start: String, rest: String) => start + rest
    }
    /** Builds a [[VcfContigHeader]] after parsing comma-delimited key-value pairs. */
    def build(implicit ctx: P[_]): P[VcfContigHeader] = csvKeyValues.flatMap { attrs: Seq[(SI, SI)] =>
      Lookup.parse(attrs).flatMap { lookup =>
        (
          lookup.andParse(key=ID, parser=name(_)) ~/
            lookup.getAndMap("length", f=_.toInt) ~/
            lookup.get("assembly")
          ).map { case (name, length, assembly) =>
          VcfContigHeader(index=0, name=name, length=length, assembly=assembly)
        }
      }
    }
    /** Parses the remaining input after the leading header key. */
    def parse[_: P]: P[VcfContigHeader] = "=<" ~/ build ~/ ">" ~ End
  }

  /** Parses the [[VcfFilterHeader]] entry/line. */
  object VcfFilterHeaderParser {
    /** Builds a [[VcfFilterHeader]] after parsing comma-delimited key-value pairs. */
    def build(implicit ctx: P[_]): P[VcfFilterHeader] = csvKeyValues.flatMap { attrs: Seq[(SI, SI)] =>
      Lookup.parse(attrs).flatMap { lookup =>
        (lookup(ID) ~/ lookup(Description)).map { case (id, description) =>
          VcfFilterHeader(id=id, description=description)
        }
      }
    }
    /** Parses the remaining input after the leading header key. */
    def parse[_: P]: P[VcfFilterHeader] = "=<" ~/ build ~/ ">" ~ End
  }

  /** Base trait for parsers for [[VcfCountAndKindHeader]]s. */
  sealed trait VcfCountAndKindHeaderParser[T<:VcfCountAndKindHeader] {
    /** Zero or more reserved header entries/lines of this type against which kind and count are validated. */
    protected def reservedHeaders: Seq[T]
    /** Mapping from reserved header entry/line ID and entry/line. */
    private val reservedHeadersMap: Map[String, T] = reservedHeaders.map { header => header.id -> header }.toMap
    /** Returns a parser for this entry/line given a lookup of attributes, id, count, and kind. */
    protected def toParser[_: P](lookup: Lookup, id: String, count: Count, kind: Kind): P[T]
    /** Parses a [[Count]] from the input for a [[VcfInfoHeader]] */
    protected def numberToCount(implicit ctx: P[_]): P[Count]
    /** Builds a header entry of type [[T]] after parsing comma-delimited key-value pairs. */
    def build[_: P]: P[T] = csvKeyValues.flatMap { attrs: Seq[(SI, SI)] =>
      Lookup.parse(attrs).flatMap { lookup =>
        // Check required fields
        (lookup(ID) ~/
          lookup.andParse(key = Number, parser = numberToCount(_)) ~/
          lookup.andParse(key = Type, parser = typeToFieldType(_))
          ).flatMap { case (id, count, kind) =>
          // Check against pre-defined/reserved INFO values
          reservedHeadersMap.get(id) match {
            case Some(info) if info.count != count =>
              lookup.si(key = Number).fail(s"expected `${info.count}` for $Number, found $count")
            case Some(info) if info.kind != kind =>
              lookup.si(key = Type).fail(s"expected `${info.kind}` for $Type, found $kind")
            case _ => toParser(lookup=lookup, id=id, count=count, kind=kind)
          }
        }
      }
    }
    /** Parses the remaining input after the leading header key. */
    def parse[_: P]: P[T] = "=<" ~/ build ~ ">" ~ End
  }

  /** Parses the [[VcfInfoHeader]] entry/line. */
  object VcfInfoHeaderParser extends VcfCountAndKindHeaderParser[VcfInfoHeader] {
    /** The reserved INFO header lines. */
    protected def reservedHeaders: Seq[VcfInfoHeader] = IndexedSeq(
      VcfInfoHeader(id="AA",        count=1,               kind=Kind.String,  description=""),
      VcfInfoHeader(id="AC",        count=OnePerAltAllele, kind=Kind.Integer, description=""),
      VcfInfoHeader(id="AD",        count=OnePerAllele,    kind=Kind.Integer, description=""),
      VcfInfoHeader(id="ADF",       count=OnePerAllele,    kind=Kind.Integer, description=""),
      VcfInfoHeader(id="ADR",       count=OnePerAllele,    kind=Kind.Integer, description=""),
      VcfInfoHeader(id="AF",        count=OnePerAltAllele, kind=Kind.Float,   description=""),
      VcfInfoHeader(id="AN",        count=1,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="BQ",        count=1,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="CIGAR",     count=OnePerAltAllele, kind=Kind.String,  description=""),
      VcfInfoHeader(id="DB",        count=0,               kind=Kind.Flag,    description=""),
      VcfInfoHeader(id="DP",        count=1,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="END",       count=1,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="H2",        count=0,               kind=Kind.Flag,    description=""),
      VcfInfoHeader(id="H3",        count=0,               kind=Kind.Flag,    description=""),
      VcfInfoHeader(id="MQ",        count=1,               kind=Kind.Float,   description=""),
      VcfInfoHeader(id="MQ0",       count=1,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="NS",        count=1,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="SB",        count=4,               kind=Kind.Integer, description=""),
      VcfInfoHeader(id="SOMATIC",   count=0,               kind=Kind.Flag,    description=""),
      VcfInfoHeader(id="VALIDATED", count=0,               kind=Kind.Flag,    description=""),
      VcfInfoHeader(id="1000G",     count=0,               kind=Kind.Flag,    description=""),
    )

    /** Parses a [[Count]] from the input for a [[VcfInfoHeader]] */
    protected def numberToCount(implicit ctx: P[_]): P[Count] = {
      Unknown.number.!.map(_ => Unknown) |
        OnePerAltAllele.number.!.map(_ => OnePerAltAllele) |
        OnePerAllele.number.!.map(_ => OnePerAllele) |
        OnePerGenotype.number.!.map(_ => OnePerGenotype) |
        CharsWhileIn("0-9").rep(1).!.map { value => Fixed(value.toInt) }
    }
    /** Returns a parser for this entry/line given a lookup of attributes, id, count, and kind. */
    protected def toParser[_: P](lookup: Lookup, id: String, count: Count, kind: Kind): P[VcfInfoHeader] = {
      (lookup(Description) ~/ lookup.get(Source) ~/ lookup.get(Version)).map {
        case (description, source, version) =>
          VcfInfoHeader(
            id          = id,
            count       = count,
            kind        = kind,
            description = description,
            source      = source,
            version     = version
          )
      }
    }
  }

  /** Parses the [[VcfFormatHeader]] entry/line. */
  object VcfFormatHeaderParser extends VcfCountAndKindHeaderParser[VcfFormatHeader] {
    /** The reserved FORMAT header lines. */
    protected def reservedHeaders: Seq[VcfFormatHeader] = IndexedSeq(
      VcfFormatHeader(id="AD",  count=OnePerAllele,    kind=Kind.Integer, description=""),
      VcfFormatHeader(id="ADF", count=OnePerAllele,    kind=Kind.Integer, description=""),
      VcfFormatHeader(id="ADR", count=OnePerAllele,    kind=Kind.Integer, description=""),
      VcfFormatHeader(id="DP",  count=1,               kind=Kind.Integer, description=""),
      VcfFormatHeader(id="EC",  count=OnePerAltAllele, kind=Kind.Integer, description=""),
      VcfFormatHeader(id="FT",  count=1,               kind=Kind.String,  description=""),
      VcfFormatHeader(id="GL",  count=OnePerGenotype,  kind=Kind.Float,   description=""),
      VcfFormatHeader(id="GP",  count=OnePerGenotype,  kind=Kind.Float,   description=""),
      VcfFormatHeader(id="GQ",  count=1,               kind=Kind.Integer, description=""),
      VcfFormatHeader(id="GT",  count=1,               kind=Kind.String,  description=""),
      VcfFormatHeader(id="HQ",  count=2,               kind=Kind.Integer, description=""),
      VcfFormatHeader(id="MQ",  count=1,               kind=Kind.Integer, description=""),
      VcfFormatHeader(id="PL",  count=OnePerGenotype,  kind=Kind.Integer, description=""),
      VcfFormatHeader(id="PP",  count=OnePerGenotype,  kind=Kind.Integer, description=""),
      VcfFormatHeader(id="PQ",  count=1,               kind=Kind.Integer, description=""),
      VcfFormatHeader(id="PS",  count=1,               kind=Kind.Integer, description=""),
    )
    /** Parses a [[Count]] from the input for a [[VcfFormatHeader]] */
    protected def numberToCount(implicit ctx: P[_]): P[Count] = {
      Unknown.number.!.map(_ => Unknown) |
        OnePerAltAllele.number.!.map(_ => OnePerAltAllele) |
        OnePerAllele.number.!.map(_ => OnePerAllele) |
        OnePerGenotype.number.!.map(_ => OnePerGenotype) |
        CharsWhileIn("0-9").rep(1).!.map { value => Fixed(value.toInt) }
    }
    /** Returns a parser for this entry/line given a lookup of attributes, id, count, and kind. */
    protected def toParser[_: P](lookup: Lookup, id: String, count: Count, kind: Kind): P[VcfFormatHeader] = {
      lookup(Description).map { description =>
        VcfFormatHeader(
          id          = id,
          count       = count,
          kind        = kind,
          description = description
        )
      }
    }
  }

  /** Parses the [[VcfAlternateAlleleHeader]] entry/line. */
  object VcfAlternateAlleleHeaderParser {
    /** Builds a [[VcfAlternateAlleleHeader]] after parsing comma-delimited key-value pairs. */
    def build(implicit ctx: P[_]): P[VcfAlternateAlleleHeader] = csvKeyValues.flatMap { attrs: Seq[(SI, SI)] =>
      Lookup.parse(attrs).flatMap { lookup =>
        (lookup(ID) ~/ lookup(Description)).map { case (id, description) =>
          VcfAlternateAlleleHeader(allele=VcfSymbolicAllele.build(id=id), description=description)
        }
      }
    }
    /** Parses the remaining input after the leading header key. */
    def parse[_: P]: P[VcfAlternateAlleleHeader] = "=<" ~/ build ~/ ">" ~ End
  }

  /** Parses the [[VcfGeneralHeader]] entry/line. */
  object VcfGeneralHeaderParser {
    def buildKeyVal[_: P](headerType: String): P[VcfGeneralHeader] = {
      CharsWhile(_ => true).!.map { value =>
        VcfGeneralHeader(headerType=headerType, id=value)
      }
    }
    /** Builds a [[VcfFilterHeader]] after parsing comma-delimited key-value pairs. */
    def buildAttr[_: P](headerType: String): P[VcfGeneralHeader] = "<" ~/ csvKeyValues.flatMap { attrs: Seq[(SI, SI)] =>
      Lookup.parse(attrs).flatMap { lookup =>
        lookup(ID).map { id =>
          val data = lookup
            .iterator
            .filter(_._1 == ID.entryName)
            .map { case (key, result) => key -> result.value }
            .toMap
          VcfGeneralHeader(
            headerType = headerType,
            id         = id,
            data       = data
          )
        }
      }
    } ~ ">"
    /** Parses the remaining input after the leading header key. */
    def parse[_: P](headerType: String): P[VcfGeneralHeader] = "=" ~/ (buildAttr(headerType) | buildKeyVal(headerType)) ~ End
  }

  /** Converts the given count to a [[Count.Fixed]] value. */
  private implicit def countToFixedVcfCount(count: Int): Count.Fixed = Count.Fixed(count)

  /** The names for the valid [[Kind]]s. */
  private val VcfFieldTypeNames: Seq[String] = Kind.values.map(_.toString)

  /** The names for the valid [[Kind]]s, not including the [[Kind.Flag]] */
  private val VcfFieldTypeNamesNoFlag: Seq[String] = {
    Kind.values.filterNot(_ == Kind.Flag).map(_.toString)
  }

  /** Parses the [[Kind]]. */
  def typeToFieldType(implicit ctx: P[_]): P[Kind] = {
    P(VcfFieldTypeNames.foldLeft[P[Unit]](Fail)(_ | _)).!.map { value =>
      Kind(value)
    }
  }

  /** Parses the [[Kind]], not including the [[Kind.Flag]] */
  def typeToFieldTypeNoFlag(implicit ctx: P[_]): P[Kind] = {
    P(VcfFieldTypeNamesNoFlag.foldLeft[P[Unit]](Fail)(_ | _)).!.map { value => Kind(value) }
  }

  // TODO: add support in fgbio for
  // - file format header line
  // - additional key-values in the various contig/INFO/FORMAT/generic header lines
  // - PEDIGREE lines? then add validate
}
