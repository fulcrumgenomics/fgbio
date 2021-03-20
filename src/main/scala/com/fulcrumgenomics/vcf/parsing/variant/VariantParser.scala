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

package com.fulcrumgenomics.vcf.parsing.variant

import com.fulcrumgenomics.fasta.SequenceMetadata
import com.fulcrumgenomics.vcf.api.Allele.{NoCallAllele, SpannedAllele, SymbolicAllele}
import com.fulcrumgenomics.vcf.api._
import com.fulcrumgenomics.vcf.parsing.util.ParseResult._
import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex
import com.fulcrumgenomics.vcf.parsing.util.ValueAndIndex.{CaptureResultAndIndex, _}
import fastparse.NoWhitespace._
import fastparse.{CharIn, CharPred, Index, P, _}

import scala.annotation.tailrec
import scala.collection.immutable.ListMap


// TODO:
// - id must be unique across variants
// - coordinate sort order...
//

sealed trait VariantParserUtil {
  def field[_ :P]: P[SI] = (Index ~ CharPred(_ != '\t').rep(1).! ~ Index).si
  def tabThenField[_ :P]: P[SI] = "\t" ~/ field
}

/** Trait to parse the contig name and return the [[SequenceMetadata]] */
sealed trait SequenceMetadataParser extends VariantParserUtil {
  def header: VcfHeader

  /** Parses the contig name and returns the associated [[SequenceMetadata]], or fails if not found. */
  def parseSequenceMetadata[_: P]: P[SequenceMetadata] = (Index ~ field ~ Index).si.flatMap { name =>
    header.dict.get(name.value) match {
      case None       => name.fail(s"contig not found in header: `${name.value}`")
      case Some(data) => success(data)
    }
  }
}

/** Trait to parse the variant position. */
sealed trait PositionParser {
  def parsePosition[_: P]: P[Int] = {
    (CharIn("1-9").! ~ CharsWhileIn("0-9").!).map { case (left, right) => (left + right).toInt }
  }
}

/** Trait to parse the variant identifier. */
sealed trait IdentifierParser {

  private def parseId[_: P]: P[SI] = {
    (Index ~ CharsWhile(char => char != ';' && char != ' ').! ~ Index).si
  }

  def parseIdentifier[_: P]: P[Seq[SI]] = {
    (Index ~ ".".! ~ Index).si.map(si => Seq(si)) | parseId.rep(1)
  }
}

/** Trait to parse the reference allele. */
sealed trait RefAlleleParser extends VariantParserUtil {

  private val DnaBases: String = "ACGTNacgtn"

  def parseRefAllele[_: P](metadata: SequenceMetadata, pos: Int): P[String] = field.flatMap { refSI =>
    val ref: String = refSI.value
    val refEnd: Int = pos + ref.length
    if (refEnd > metadata.length) {
       refSI.fail(s"variant end was past the end of the reference: `$refEnd > ${metadata.length}`")
    }
    else {
      ref.indexWhere(base => !DnaBases.contains(base)) match {
        case -1    => success(ref)
        case index =>
          val prefix      = ref.take(index)
          val invalidBase = ref(index)
          val suffix      = ref.drop(index + 1)
          refSI.fail(s"invalidate ref base: `$prefix[$invalidBase]$suffix")
      }
    }
  }
}

/** Trait to parse the alternate allele(s). */
sealed trait AltAlleleParser extends VariantParserUtil {

  def dna[_: P]: P[Allele] = CharIn("ACGTNacgtn").rep(1).!.map(Allele.apply)
  def spanendAllele[_: P]: P[Allele] = "*".!.map(_ => SpannedAllele)
  def missingAllele[_: P]: P[Allele] = ".".!.map(_ => NoCallAllele)
  def symbolicAllele[_: P]: P[Allele] = ("<" ~ CharPred(_ != '>').! ~ ">").map(SymbolicAllele.apply)
  // TODO: breakend

  def alt[_: P]: P[Allele] = dna | spanendAllele | missingAllele | symbolicAllele

  def alts[_: P]: P[Seq[Allele]] = alt.rep(min=1, sep=",")

  def parseAltAlleles[_: P](ref: String): P[AlleleSet] = alts.map { alleles =>
    AlleleSet(
      ref  = Allele(ref),
      alts = alleles
    )
  }
}

/** Trait to parse the variant quality. */
sealed trait QualParser {
  private def double[_: P]: P[Double] = {
    (CharIn("1-9") ~ CharsWhileIn("0-9").rep ~ ("." ~ CharsWhileIn("0-9")).?.!).!.map(_.toDouble)
  }
  def parseQual[_: P]: P[Option[Double]] = ".".!.map(_ => None) | double.map(Some(_))
}

/** Trait to parse the variant filter(s). */
sealed trait FilterParser {
  def header: VcfHeader

  private val illegalFilterChars: String = "; \t\n"

  private def value[_: P]: P[SI] = (Index ~ CharPred(!illegalFilterChars.contains(_)).rep(1).! ~ Index).si

  private def filterValue[_: P]: P[VcfFilterHeader] = {
    value.flatMap { key =>
      header.filter.get(key.value) match {
        case Some(filter) => success(filter)
        case None         => key.fail(s"FILTER was not found in the header: `${key.value}`")
      }
    }
  }
  private def filterValues[_: P]: P[Option[Seq[VcfFilterHeader]]] = filterValue.rep(min=0, sep=":").map(Some(_))
  private def missingValue[_: P]: P[Option[Seq[VcfFilterHeader]]] = ".".!.map(_ => None)
  private def pass[_: P]: P[Option[Seq[VcfFilterHeader]]] = "PASS".!.map(_ => Some(Seq.empty))
  def parseFilter[_: P]: P[Option[Seq[VcfFilterHeader]]] = pass | missingValue | filterValues
}

sealed trait FieldUtil {
  protected case class Field[T](key: SI, field: VcfCountAndKindHeader, values: Seq[T] = Seq.empty)

  private sealed trait ParseResult {
    def isSuccess: Boolean
    def isFailure: Boolean = !isSuccess
    def success: ParseSuccess
    def failure: ParseFailure
  }
  private case class ParseSuccess(headerField: VcfCountAndKindHeader, results: Seq[Any]) extends ParseResult {
    override def isSuccess: Boolean = true
    def success: ParseSuccess = this
    def failure: ParseFailure = throw new NoSuchElementException
  }
  private case class ParseFailure(headerField: VcfCountAndKindHeader, value: SI, ex: Exception) extends ParseResult {
    override def isSuccess: Boolean = false
    def success: ParseSuccess = throw new NoSuchElementException
    def failure: ParseFailure = this
  }

  @tailrec
  private def parseRecursive(values: Seq[SI], curSuccess: ParseSuccess): ParseResult = {
    values match {
      case Seq()                  => curSuccess
      case Seq(value, tail @ _*) =>
        val parseResult = try {
          curSuccess.copy(
            results = curSuccess.results :+ curSuccess.headerField.kind.parse(value.value)
          )
        } catch {
          case ex: Exception => ParseFailure(headerField=curSuccess.headerField, value=value, ex=ex)
        }
        // Developer note: this is to make the function tail recursive
        parseResult match {
          case failure: ParseFailure => failure
          case success: ParseSuccess => parseRecursive(values=tail, curSuccess=success)
        }
    }
  }

  private def parseValues(headerField: VcfCountAndKindHeader, values: Seq[SI]): ParseResult = {
    val startSuccess = ParseSuccess(headerField=headerField, results=Seq.empty)
    values match {
      case Seq(value) if value.value == "." =>
        parseRecursive(values=Seq.empty, curSuccess=startSuccess)
      case _         =>
        parseRecursive(values=values, curSuccess=startSuccess)
    }
  }

  @tailrec
  private def parseKindRecursive(inputFields: Seq[Field[SI]], pastResults: Seq[ParseResult]): Seq[ParseResult] = {
    inputFields match {
      case Seq()                => pastResults
      case Seq(data, tailFields @ _*) =>
        parseValues(headerField=data.field, values=data.values) match {
          case failure: ParseFailure => Seq(failure)
          case success: ParseSuccess => parseKindRecursive(inputFields=tailFields, pastResults=pastResults :+ success)
        }
    }
  }

  protected def parseKind[_: P](inputFields: Seq[Field[SI]], fieldName: String): P[Seq[(String,Any)]] = {
    parseKindRecursive(inputFields=inputFields, pastResults=Seq.empty) match {
      case Seq(ParseFailure(headerField, value, ex)) =>
        value.fail(s"could not parse $fieldName `${value.value}` as type `${headerField.kind}`: ${ex.getMessage}")
      case results =>
        success(results.map(_.success).map {
          case ParseSuccess(headerField, Seq(result)) => headerField.id -> result
          case ParseSuccess(headerField, results)     => headerField.id -> results
        })
    }
  }

  protected def parseCount[_: P](inputFields: Seq[Field[SI]], alleles: AlleleSet, fieldName: String): P[Seq[Field[SI]]] = {
    inputFields.find { case Field(_, field, values) =>
      val countOption: Option[Int] = field.count match {
        case VcfCount.Fixed(count)    => Some(count)
        case VcfCount.OnePerAllele    => Some(alleles.size)
        case VcfCount.OnePerAltAllele => Some(alleles.size - 1)
        case VcfCount.Unknown         => None
        case VcfCount.OnePerGenotype  => throw new Exception(s"bug: found one-per-genotype in $fieldName.${field.id}")
      }
      countOption.exists(count => count != values.length)
    } match {
      case None => success(inputFields)
      case Some(Field(key, field, values)) =>
        values.headOption.getOrElse(key).fail(
          s"$fieldName key `${field.id}`: expected ${field.count} values but found ${values.length} values"
        )
    }
  }

  protected def parseDuplicateKeys[_: P](inputFields: Seq[Field[SI]], fieldName: String): P[Seq[Field[SI]]] = {
    inputFields.groupBy(_.key.value).find(_._2.length > 1) match {
      case None                => success(inputFields)
      case Some((key, fields)) => fields.head.key.fail(f"found $fieldName key `$key` ${fields.length} times")
    }
  }
}

sealed trait InfoParser extends FieldUtil {
  def header: VcfHeader

  private def key[_: P]: P[SI] = (Index ~ ("1000G".! | (CharIn("A-Za-z_") ~ CharIn("0-9A-Za-z_.").rep).!) ~ Index).si

  /** Parses the INFO key, making sure it is found in the VCF header.  Must be `^([A-Za-z ][0-9A-Za-z .]*|1000G)$,`. */
  def infoKey[_: P]: P[Field[SI]] = {
    key.flatMap { key =>
      this.header.info.get(key.value) match {
        case None       => key.fail(f"INFO key not found in the header: `${key.value}`")
        case Some(info) => success(Field[SI](key, info))
      }
    }
  }

  /** Parses the value for an INFO field. */
  def infoValue[_: P]: P[SI] = (Index ~ CharPred(c => c != ';' && c != '=' && c != ',' && c != '\t').rep(1).! ~ Index).si

  /** Parses an INFO key and value (may be a list) */
  def infoKeyAndValues[_: P]: P[Field[SI]] = {
    (infoKey ~ ("=" ~ infoValue.rep(min=1, sep=",")).?).map {
      case (data: Field[SI], None) => data
      case (data, Some(values))    => data.copy(values=values)
    }
  }

  def infos[_: P]: P[Seq[Field[SI]]] = ".".!.map(_ => Seq.empty) | infoKeyAndValues.rep(sep=";")

  private def parseMissing[_: P](missing: Field[SI], infoValues: Seq[Field[SI]]): P[ListMap[String, Any]] = {
    if (infoValues.length != 1) {
      missing.key.fail("found `.` with other INFO fields")
    }
    else if (missing.values.nonEmpty) {
      missing.values.head.fail("INFO field `.` cannot have a value")
    }
    else {
      success(ListMap.empty)
    }
  }

  def parseInfos[_: P](alleles: AlleleSet): P[ListMap[String, Any]] = infos.flatMap { infoValues=>
    infoValues.find(_.key.value == ".") match {
      case Some(missing) => parseMissing(missing=missing, infoValues=infoValues)  // empty INFO
      case None          => // non-empty INFOs
        parseDuplicateKeys(inputFields=infoValues, fieldName="INFO")
          .flatMap { _ => parseCount(inputFields=infoValues, alleles=alleles, fieldName="INFO") }
          .flatMap { _ => parseKind(inputFields=infoValues, fieldName="INFO") }
          .map { values =>
            val builder = ListMap.newBuilder[String, Any]
            builder ++= values
            builder.result()
          }
    }
  }
}

sealed trait FormatKeysParser {
  def header: VcfHeader

  private def checkKeyInHeader[_: P](key: SI): P[VcfFormatHeader] = {
    this.header.format.get(key.value) match {
      case None         => key.fail(f"FORMAT key not found in the header: `${key.value}`")
      case Some(format) => success(format)
    }
  }

  def gtKey[_: P]: P[VcfFormatHeader] = (Index ~ "GT".! ~ Index).si.flatMap {
    checkKeyInHeader(_)
  }

  def rawKey[_: P]: P[SI] = (Index ~ (CharIn("A-Za-z_") ~ CharIn("0-9A-Za-z_.").rep).! ~ Index).si

  def rawKeys[_: P]: P[Seq[SI]] =  ".".!.map(_ => Seq.empty) | rawKey.rep(1, sep=":")

  def parseFormatKeys[_: P]: P[Seq[VcfFormatHeader]] = {
    rawKeys.flatMap { keys: Seq[SI] =>
      if (keys.isEmpty) success(keys)
      else if (keys.head.value != "GT") keys.head.fail(s"the first FORMAT key must be GT")
      else { // Check for duplicate keys
        keys.groupBy(_.value).find(_._2.length > 1) match {
          case None => success(keys)
          case Some((key, values)) => values.head.fail(f"found FORMAT key `$key` ${values.length} times")
        }
      }
    }.flatMap { keys =>
      // Find each key in the header
      keys.find(key => !this.header.format.contains(key.value)) match {
        case None      => success(keys.map(key => this.header.format(key.value)))
        case Some(key) => key.fail(f"FORMAT key not found in the header: `${key.value}`")
      }
    }
  }
}

sealed trait GenotypeParser extends VariantParserUtil with FieldUtil {
  def header: VcfHeader

  def gtValue[_: P]: P[SI] = (Index ~ ("0".! | (CharIn("1-9") ~ CharIn("0-9").rep).!) ~ Index).si

  def gt[_: P](alleles: AlleleSet): P[Allele] = ".".!.map(_ => NoCallAllele) | gtValue.flatMap { v: SI =>
    val index = v.value.toInt
    alleles.get(index) match {
      case None         => v.fail(s"genotype index `$index` is out of range: [0, ${alleles.size}]")
      case Some(allele) => success(allele)
    }
  }
  def parseGT[_: P](alleles: AlleleSet): P[(Seq[Allele], Boolean)] = {
    gt(alleles).rep(min=1, sep="/").map(alleles => (alleles, true)) |
      gt(alleles).rep(min=1, sep="|").map(alleles => (alleles, false))
  }

  /** Parses the value for an FORMAT field. */
  def genotypeValue[_: P]: P[SI] = (Index ~ CharPred(c => c != ':' && c != ',' && c != '\t').rep(1).! ~ Index).si

  /** Parses an FORMAT key and value (may be a list) */
  def genotypeValues[_: P]: P[Seq[SI]] = genotypeValue.rep(min=0, sep=",")

  def genotype[_: P](formatKeys: Seq[VcfFormatHeader], alleles: AlleleSet): P[Genotype] = {
    require(formatKeys.headOption.forall(_.id == "GT"))
    val tailFormatKeys = formatKeys.tail
    (".".!.map(_ => (Seq.empty, false)) | parseGT(alleles=alleles)).flatMap { case (calls: Seq[Allele], phased: Boolean) =>
      ":" ~ genotypeValues.rep(min=0, max=formatKeys.length-1, sep=":").flatMap { fieldValues: Seq[Seq[SI]] =>
        if (fieldValues.length > tailFormatKeys.length) {
          fieldValues(tailFormatKeys.length).head.fail(f"found too many FORMAT values ${fieldValues.length} > ${tailFormatKeys.length}")
        }
        else { // Build some fields for parsing
          // Unspecified trailing fields are allowed
//          val missingFields: Seq[Field[SI]] = tailFormatKeys.drop(fieldValues.length).map { field =>
//            val key: SI = ValueAndIndex(value=field.id, start=0, end=0) // FIXME: start and end
//            Field(key=key, field=field, values=Seq.empty) // FIXME: what should these be?
//          }
          // Specified fields
          val fields: Seq[Field[SI]] = fieldValues.zip(tailFormatKeys).map { case (attrValues: Seq[SI], field: VcfFormatHeader) =>
            val key: SI = ValueAndIndex(value=field.id, start=0, end=0) // FIXME: start and end
            Field(key=key, field=field, values=attrValues)
          }
//          assert(missingFields.length + fields.length == tailFormatKeys.length)
//          success(fields ++ missingFields)
          success(fields)
        }
      }
      .flatMap { fields: Seq[Field[SI]] => parseCount(inputFields=fields, alleles=alleles, fieldName="FORMAT") }
      .flatMap { fields: Seq[Field[SI]] => parseKind(inputFields=fields, fieldName="FORMAT") }
      .map { formatValues: Seq[(String, Any)] =>
        val builder = ListMap.newBuilder[String, Any]
        builder ++= formatValues
        builder.result()
      }.map { attributes: Map[String, Any] =>
        Genotype(
          alleles = alleles,
          sample  = "", // NB: this will be updated later
          calls   = calls.toIndexedSeq,
          phased  = phased,
          attrs   = attributes
        )
      }
    }
  }

  def missingGenotype[_: P](alleles: AlleleSet): P[Genotype] = {
    ".".!.map(_ => Genotype(alleles=alleles, sample="", calls=IndexedSeq.empty))
  }

  def parseGenotype[_: P](formatKeys: Seq[VcfFormatHeader], alleles: AlleleSet): P[Genotype] = {
    (missingGenotype(alleles=alleles) | genotype(formatKeys=formatKeys, alleles=alleles))
  }
}

class VariantParser(val header: VcfHeader)
  extends SequenceMetadataParser
    with PositionParser
    with IdentifierParser
    with RefAlleleParser
    with AltAlleleParser
    with QualParser
    with FilterParser
    with InfoParser
    with FormatKeysParser
    with GenotypeParser {

  private val numSamples: Int = header.samples.length

  private def parseChromPosAndIdentifier[_: P]: P[(SequenceMetadata, Int, Seq[SI])] = {
    Start ~/
      parseSequenceMetadata ~/ // CHROM
      "\t" ~/ parsePosition ~/ // POS
      "\t" ~/ parseIdentifier ~/  // ID
      "\t"
  }

  def variant[_: P]: P[Variant] = {
    parseChromPosAndIdentifier.flatMap { case (metadata, pos, id) =>
      parseRefAllele(metadata=metadata, pos=pos).map { ref => (metadata, pos, id, ref) } ~/ "\t"
    }.flatMap { case (metadata, pos, id, ref) =>
      (
        parseAltAlleles(ref=ref) ~/ // ALT
        "\t" ~/ parseQual ~/        // QUAL
        "\t" ~/ parseFilter ~/      // FILTER
        "\t"
      ).map { case (alleleSet, qual, filters) =>
        (metadata, pos, id, ref, alleleSet, qual, filters)
      }
    }.flatMap { case (metadata, pos, id, ref, alleleSet, qual, filters) =>
      // INFO
      parseInfos(alleleSet).map { infos =>
        (metadata, pos, id, ref, alleleSet, qual, filters, infos)
      }
    }.flatMap { case (metadata, pos, id, ref, alleleSet, qual, filters, infos) =>
      if (this.numSamples == 0) {
        "" ~ End.map { _ =>
          (metadata, pos, id, ref, alleleSet, qual, filters, infos, Seq.empty[Genotype])
        }
      }
      else {
        ("\t" ~/ parseFormatKeys).flatMap { formats: Seq[VcfFormatHeader] =>
          "\t" ~/
            parseGenotype(formatKeys=formats, alleles=alleleSet)
              .rep(exactly=this.numSamples, sep="\t")
              .map { genotypes => (metadata, pos, id, ref, alleleSet, qual, filters, infos, genotypes) } ~ End
        }
      }
    }.map { case (metadata, pos, id, _, alleleSet, qual, filters, infos, genotypes) =>
      val filterStrings = filters match {
        case None           => Seq()
        case Some(Seq())    => Seq("PASS")
        case Some(_filters) => _filters.map(_.id)
      }

      new Variant(
        chrom     = metadata.name,
        pos       = pos,
        id        = id.map(_.value),
        alleles   = alleleSet,
        qual      = qual,
        filters   = filterStrings.toSet,
        attrs     = infos,
        genotypes = header.samples.zip(genotypes).toMap
      )
    }
  }
}