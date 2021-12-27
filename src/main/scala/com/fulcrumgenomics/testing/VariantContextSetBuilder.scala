/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.testing

import java.nio.file.Files
import java.util.Collections

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.fasta.SequenceDictionary
import htsjdk.variant.variantcontext._
import htsjdk.variant.variantcontext.writer.{Options, VariantContextWriterBuilder}
import htsjdk.variant.vcf.{VCFFileReader, VCFHeader, VCFHeaderLine}

import scala.collection.JavaConverters._
import scala.collection.compat._
import scala.collection.mutable.ListBuffer
import scala.reflect.ClassTag

object VariantContextSetBuilder {
  @deprecated(since="0.9.0", message="Use com.fulcrumgenomics.testing.VcfBuilder instead.")
  def apply(sampleName: String): VariantContextSetBuilder = {
    new VariantContextSetBuilder(List(sampleName))
  }
}

/**
  * A class to make artificial variant contexts and VCFs just a little bit easier to generate.
  *
  * This builder uses the default sequence dictionary from [SamBuilder] by default.
  */
@deprecated(since="0.9.0", message="Use com.fulcrumgenomics.testing.VcfBuilder instead.")
class VariantContextSetBuilder(sampleNames: Seq[String] = List("Sample")) extends Iterable[VariantContext] {

  if (sampleNames.isEmpty) throw new IllegalArgumentException("At least one sample name must be given")

  private var _header = {
    val h = new VCFHeader(Collections.emptySet[VCFHeaderLine](), sampleNames.toList.asJava)
    h.setSequenceDictionary(new SamBuilder().header.getSequenceDictionary)
    h
  }

  private val variants = new ListBuffer[VariantContext]()

  def header: VCFHeader = this._header

  lazy val dict: SequenceDictionary = {
    import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
    this.header.getSequenceDictionary.fromSam
  }

  /** Sets the VCF header for this builder.  This should be set before adding variants. */
  def setHeader(header: VCFHeader): this.type = yieldingThis(this._header = header)

  /** Sets the sequence dictionary for this builder.  This should be set before adding variants. */
  def setSequenceDictionary(dict: SequenceDictionary): this.type = {
    import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
    yieldingThis(this._header.setSequenceDictionary(dict.asSam))
  }

  /** Adds the header line to the header for this builder.  This should be set before adding variants. */
  def addMetaDataLine(headerLine: VCFHeaderLine): this.type = yieldingThis(this._header.addMetaDataLine(headerLine))

  /** Adds a variant to the builder.
    *
    * If a variant context already exists with the same contig, start, stop, and variant alleles, then we add this
    * genotype to that variant, otherwise we add a new variant and genotype.
    *
    * The variant alleles must always be given, as we use them to determine if this variant already exists when adding
    * a genotype for a new sample.
    *
    * If no genotype alleles are given, the genotype is a no call.
    *
    * If the sample name is not given, the first sample name in the header (in order) is used.  There must not already
    * be a genotype for this variant with the given sample name, and the sample name must be present in the header.
    *
    * Genotype attributes may be given with the `genotypeAttributes` parameter.  The value for the `GQ` and `DP`
    * attributes must have type [[Int]].  The value for the `AD` and `PL` attributes must have either type
    * [[IterableOnce]] or [[Array]], with each item having type [[Int]].  For example:
    * {{{
    *   val builder = new VariantContextSetBuilder()
    *   builder.addVariant(
    *     start=1,
    *     variantAlleles=Seq("A", "C"),
    *     genotypeAlleles=Seq("A", "C"),
    *     genotypeAttributes=Map("DP" -> 10, "AD" -> Array(6, 4))
    *   )
    * }}}
    *
    * */
  def addVariant(refIdx: Int = 0,
                 start: Long,
                 variantAlleles: List[String],
                 variantAttributes: Map[String,Any] = Map.empty,
                 genotypeAlleles: List[String] = List.empty,
                 genotypeAttributes: Map[String,Any] = Map.empty,
                 sampleName: Option[String] = None,
                 phased: Boolean = false): this.type = {
    if (!sampleName.forall { sn => this.header.getGenotypeSamples.contains(sn)}) {
      throw new IllegalArgumentException(s"Sample with name $sampleName not found in the VCF header.")
    }
    if (variantAlleles.isEmpty) {
      throw new IllegalArgumentException("No alleles given")
    }
    if (!genotypeAlleles.forall(a => a == Allele.NO_CALL_STRING || variantAlleles.contains(a))) {
      throw new IllegalArgumentException("A genotype allele not found in variant alleles")
    }
    val contig  = this.dict(refIdx).name
    val alleles = toAlleles(variantAlleles)
    val stop    = VariantContextUtils.computeEndFromAlleles(alleles.asJava, start.toInt, -1)
    // check to see if there are already genotypes for this variant
    val (ctxBuilder, prevGenotypes) = this.variants.find { ctx =>
      ctx.getContig == contig &&
        ctx.getStart == start &&
        ctx.getEnd == stop &&
        ctx.getAlleles.asScala.map(_.getBaseString).sorted.mkString == variantAlleles.sorted.mkString // only an approximation
    } match {
      case Some(ctx) =>
        // remove this variant context, and create a builder based off of it, and get the current genotypes.
        this.variants -= ctx
        (new VariantContextBuilder(ctx), ctx.getGenotypesOrderedByName.toList)
      case None      =>
        // create a new variant context builder and an empty set of genotypes
        (new VariantContextBuilder("source", contig, start, stop, alleles.asJava), List.empty[Genotype])
    }
    // get the reference allele
    val referenceAllele = ctxBuilder.getAlleles.find(_.isReference)
    // create a genotype builder.
    val name = sampleName.getOrElse(header.getSampleNamesInOrder.iterator().next())
    val genotypeBuilder = genotypeAlleles match {
      case Nil      => new GenotypeBuilder(name, Collections.singletonList(Allele.create(Allele.NO_CALL_STRING)))
      case gAlleles => new GenotypeBuilder(name, toAlleles(gAlleles, referenceAllele=referenceAllele).asJava)
    }
    genotypeBuilder.phased(phased)
    // Set the genotype attributes, with GQ/DP/AD/PL being special-cased
    // NB: no check is performed that attributes that are per-allele match the given # of alleles.
    genotypeAttributes.foreach {
      case ("GQ", v) => genotypeBuilder.GQ(anyToType[Integer]("GQ", v, "Int"))
      case ("DP", v) => genotypeBuilder.DP(anyToType[Integer]("DP", v, "Int"))
      case ("AD", v) => genotypeBuilder.AD(anyToArray[Integer]("AD", v, "Int").map(_.toInt))
      case ("PL", v) => genotypeBuilder.PL(anyToArray[Integer]("PL", v, "Int").map(_.toInt)) // NB: we don't support Doubles
      case (k,v)     => genotypeBuilder.attribute(k, v)
    }
    val genotype = genotypeBuilder.make()
    // check the sample doesn't already exists.
    prevGenotypes.find { g => g.getSampleName == genotype.getSampleName } match {
      case Some(g) => throw new IllegalArgumentException(s"Genotype already exists for the sample ${g.getSampleName}: ${g.toBriefString}")
      case None    =>
    }
    val genotypes = genotype +: prevGenotypes
    ctxBuilder.genotypes(genotypes:_*)
    ctxBuilder.attributes(variantAttributes.asJava)
    yieldingThis(this.variants.append(ctxBuilder.make()))
  }

  /** Returns `v` as type `T`.  Throws an exception if `v` is not assignable to type `T`.*/
  private def anyToType[T : ClassTag](k: String, v: Any, typeName: String): T = {
    val clazz = implicitly[ClassTag[T]].runtimeClass
    v match {
      case value if clazz.isInstance(v) => value.asInstanceOf[T]
      case _ => throw new IllegalArgumentException(s"$k attribute value must be an $typeName, found ${v.getClass.getSimpleName}")
    }
  }

  /** Returns `v` as type `Array[T]`.  Throws an exception if `v` is not assignable to `TraversableOnce[_]` or `Array[_]`,
    * or if any element in `v` is not of type `T` .*/
  private def anyToArray[T : ClassTag](k: String, v: Any, typeName: String): Array[T] = {
    val values = v match {
      case _v: IterableOnce[_] => _v.iterator.toArray[Any]
      case _v: Array[_]           => _v
      case _                      =>
        throw new IllegalArgumentException(s"$k attribute value must be an IterableOnce[$typeName], found ${v.getClass.getSimpleName}")
    }
    values.map(anyToType[T](k, _, typeName))
  }

  /** Gets an iterator over the variants. */
  override def iterator: Iterator[VariantContext] = this.variants.iterator

  /** Returns the number of variants added to this builder. */
  override def size: Int = this.variants.size

  /** Writes the contents of the record set to the provided file path. */
  def write(path: PathToVcf) : Unit = {
    val writer = new VariantContextWriterBuilder()
      .setReferenceDictionary(this.header.getSequenceDictionary)
      .setOutputFile(path.toFile)
      .setOption(Options.INDEX_ON_THE_FLY)
      .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
      .build()
    writer.writeHeader(_header)
    this.iterator.foreach(writer.add)
    writer.close()
  }

  /** Writes the contents to a temporary file that will be deleted when the JVM exits. */
  def toTempFile(deleteOnExit: Boolean = true): PathToVcf = {
    val path = Files.createTempFile("VariantContextSet.", ".vcf.gz")
    if (deleteOnExit) path.toFile.deleteOnExit()
    write(path)
    path
  }

  /** Writes the contents to a temporary file that will be deleted when the JVM exits and
    * provides a [VCFFileReader]. */
  def toVcfFileReader(deleteOnExit: Boolean = true): VCFFileReader = {
    val path = toTempFile(deleteOnExit=deleteOnExit)
    new VCFFileReader(path.toFile, true)
  }

  /** Given a list of alleles as strings, typically for a genotype, gives a list of alleles of type [Allele].  If
    * a reference allele is given, then we will use that to set if an input allele is reference, otherwise we assume
    * it's the first allele in the input list.
    */
  private def toAlleles(alleles: List[String], referenceAllele: Option[Allele] = None): List[Allele] = {
    referenceAllele match {
      case Some(r) =>
        alleles.map { allele => Allele.create(allele, r.getBaseString == allele) }
      case None    =>
        alleles.zipWithIndex.map { case (allele: String, idx: Int) => Allele.create(allele, idx == 0) }
    }
  }

  private def yieldingThis(f: => Unit): this.type = { f; this}
}
