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

package com.fulcrumgenomics.testing

import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.vcf.api.Allele.NoCallAllele
import com.fulcrumgenomics.vcf.api._

import scala.collection.compat._
import scala.collection.immutable.ListMap
import scala.collection.mutable
import scala.util.Try

object VcfBuilder {
  /** The default header that is used when making a new [[VcfBuilder]]. */
  val DefaultHeader: VcfHeader = {
    val contigs = new SamBuilder().dict.getSequences
      .map(s => VcfContigHeader(s.getSequenceIndex, s.getSequenceName, length=Some(s.getSequenceLength), assembly=Option(s.getAssembly)))
      .toIndexedSeq

    VcfHeader(
      contigs = contigs,
      infos   = Seq(
        VcfInfoHeader(id="AC",  count=VcfCount.OnePerAltAllele, kind=VcfFieldType.Integer, description="Alternate allele counts in genotypes."),
        VcfInfoHeader(id="DP",  count=VcfCount.Fixed(1),        kind=VcfFieldType.Integer, description="Depth across all samples.")
      ),
      formats = Seq(
        VcfFormatHeader(id="GT", count=VcfCount.Fixed(1),       kind=VcfFieldType.String, description="Genotype string."),
        VcfFormatHeader(id="AD", count=VcfCount.OnePerAllele,   kind=VcfFieldType.Integer, description="Depth per allele."),
        VcfFormatHeader(id="GQ", count=VcfCount.Fixed(1),       kind=VcfFieldType.Integer, description="Genotype quality."),
        VcfFormatHeader(id="PL", count=VcfCount.OnePerGenotype, kind=VcfFieldType.Integer, description="Phred scaled genotype likelihoods.")
      ),
      filters = Seq(
        VcfFilterHeader(id="LowQD", description="Low Quality/Depth value"),
        VcfFilterHeader(id="LowAB", description="Low/poor allele balance.")
      ),
      others  = Seq(
        VcfGeneralHeader(headerType="ALT", id="NON_REF", data=Map("Description" -> "Represents any non-reference allele."))
      ),
      samples = IndexedSeq()
    )
  }

  /**
    * Convenience case class used to build up genotypes.
    *
    * @param sample the name of the sample being genotyped
    * @param gt a genotype either in VCF format (e.g. 0/1) or using actual alleles (e.g. A/T)
    * @param dp optional depth (DP)
    * @param ad optional set of allele depths (AD)
    * @param pl optional set of phred-scale genotype likelihoods (PL)
    * @param attrs extended attributes (must be in FORMAT headers)
    */
  case class Gt(sample: String,
                gt: String,
                dp: Int = -1,
                ad: Seq[Int] = Seq.empty,
                pl: Seq[Int] = Seq.empty,
                attrs: Map[String,Any] = Map.empty
               ) {

    def allAttrs: Map[String, Any] = {
      attrs ++
        (if (dp < 0) Map.empty else Map("DP" -> dp)) ++
        (if (ad.isEmpty) Map.empty else Map("AD" -> ad.toIndexedSeq)) ++
        (if (pl.isEmpty) Map.empty else Map("DP" -> pl.toIndexedSeq))
    }
  }

  /** Constructs a VcfBuilder using the supplied header. */
  def apply(header: VcfHeader): VcfBuilder = new VcfBuilder(header)

  /** Constructs a VcfBuilder using the [[DefaultHeader]] and the set of samples provided. */
  def apply(samples: Seq[String]): VcfBuilder = {
    require(samples.distinct.size == samples.size, s"${samples.mkString(",")} contains duplicate sample names.")
   this.apply(DefaultHeader.copy(samples=samples.toIndexedSeq))
  }
}

/** Class for building VCFs for testing purposes. */
class VcfBuilder private (initialHeader: VcfHeader) extends Iterable[Variant] {
  /** Genomic location as (sequence_index, position). */
  private type Location = (Int, Int)
  private var _header: VcfHeader = initialHeader
  private val variants: mutable.Map[Location, Variant] = mutable.TreeMap()

  /** Provides access to the header that will be written to the VCF. */
  def header: VcfHeader = _header

  /** Adds one or more INFO headers to the VCF Header. */
  def withInfoHeaders(headers: VcfInfoHeader*): this.type = {
    this._header = this._header.copy(infos = this._header.infos ++ headers)
    this
  }

  /** Adds one or more FORMAT headers to the VCF Header. */
  def withFormatHeaders(headers: VcfFormatHeader*): this.type = {
    this._header = this._header.copy(formats = this._header.formats ++ headers)
    this
  }

  /** Adds one or more FILTER headers to the VCF Header. */
  def withFilterHeaders(headers: VcfFilterHeader*): this.type = {
    this._header = this._header.copy(filters = this._header.filters ++ headers)
    this
  }

  /** Adds one or more non-INFO/FORMAT/FILTER headers to the VCF Header. */
  def withOtherHeaders(headers: VcfGeneralHeader*): this.type = {
    this._header = this._header.copy(others = this._header.others++ headers)
    this
  }

  /** Adds a variant to the set being built.  The variant should contain all information required as it is not
    * possible to update a variant once added. If a variant already exists at the given position an exception
    * will be thrown.
    *
    * The genotypes are specified using instance of the [[Gt]] helper class.  The genotype strings within the [[Gt]]
    * objects may be either numeric like in a VCF (e.g. `0/1`) or use allele strings (e.g. `A/T`.)
    *
    * The variant must also be valid, e.g. not reference INFO, FILTER or FORMAT fields that are no in the header,
    * and not have alleles in genotypes that are not in the set of alleles for the variant.
    *
    * If the header contains more samples than there are genotype give, simple diploid no-call genotypes will
    * be inserted for the remaining samples.
    */
  def add(chrom: String = this.header.contigs.head.name,
          pos: Int,
          id: String = ".",
          alleles: Seq[String],
          qual: Int = -1,
          info: Map[String, Any] = Map.empty,
          filters: IterableOnce[String] = Set.empty,
          gts: Seq[Gt] = Seq.empty
         ): this.type = {

    require(_header.dict.getSequenceIndex(chrom) >= 0, s"Unknown chrom: $chrom")
    val key = (_header.dict.getSequenceIndex(chrom), pos)

    // Make sure things are relatively valid
    require(!variants.contains(key), s"Variant already exists at position $chrom:$pos")
    require(alleles.nonEmpty, s"Must specify at least one allele.")
    info.keys.foreach(k => require(this._header.info.contains(k), s"No INFO header for key $k"))
    filters.foreach(f => require(this._header.filter.contains(f), s"No FILTER header for key $f"))
    require(gts.map(_.sample).toSet.size == gts.size, s"Non-unique sample names in genotypes.")

    val alleleSet = AlleleSet(alleles:_*)

    val calledGenotypes = gts.map { g =>
      val phased    = g.gt.contains("|")
      val separator = if (phased) '|' else '/'
      val calls     = g.gt.split(separator).map(s => toAllele(s, alleleSet)).toIndexedSeq

      Genotype(
        alleles = alleleSet,
        sample  = g.sample,
        calls   = calls,
        phased  = phased,
        attrs   = g.allAttrs
      )
    }

    val noCalls = _header.samples.diff(calledGenotypes.map(_.sample)).map { s =>
      Genotype(alleles=alleleSet, sample=s, calls=IndexedSeq(NoCallAllele, NoCallAllele))
    }

    variants.put(key, Variant(
      chrom     = chrom,
      pos       = pos,
      id        = if (id == ".") None else Some(id),
      alleles   = alleleSet,
      qual      = if (qual < 0) None else Some(qual),
      attrs     = ListMap(info.toSeq:_*),
      filters   = filters.toSet,
      genotypes = (calledGenotypes ++ noCalls).map(gt => gt.sample -> gt).toMap
    ))

    this
  }

  /** Converts either a numeric index or an allele string into an Allele. */
  private def toAllele(s: String, alleles: AlleleSet): Allele = {
    Try { alleles(s.toInt) }
      .recover { case _ => Allele(s) }
      .get
  }

  /** Returns an iterator over the variants in chromosomal order. */
  def iterator: Iterator[Variant] = this.variants.valuesIterator

  /** Writes the contents of the record set to the provided file path. */
  def write(path: PathToVcf) : Unit = {
    val writer = VcfWriter(path, this._header)
    writer ++= iterator
    writer.close()
  }

  /** Writes the contents to a temporary file that will be deleted when the JVM exits. */
  def toTempFile(deleteOnExit: Boolean = true): PathToVcf = {
    val path = Files.createTempFile("VcfRecordSet.", ".vcf.gz")
    if (deleteOnExit) path.toFile.deleteOnExit()
    write(path)
    path
  }

  /** Creates a VcfSource over the records stored in a temporary file. */
  def toSource: VcfSource = VcfSource(toTempFile())
}
