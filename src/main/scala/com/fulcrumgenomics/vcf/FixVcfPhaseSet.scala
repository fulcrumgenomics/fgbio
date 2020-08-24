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

package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, PathToVcf, SafelyClosable}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, ProgressLogger}
import com.fulcrumgenomics.vcf.api._
import enumeratum.EnumEntry

import scala.collection.{immutable, mutable}


@clp(group=ClpGroups.VcfOrBcf, description=
  """
    |Adds/fixes the phase set (PS) genotype file.
    |
    |The VCF specification allows phased genotypes to be annotated with the `PS` (phase set) `FORMAT` field.  The value
    |should be a non-negative integer, corresponding to the position of the first variant in the phase set.  Some tools
    |will output a non-integer value, as well as describe this field as having non-Integer type in the VCF header.  This
    |tool will update the phase set (`PS`) `FORMAT` field to be VCF spec-compliant.  Phased genotypes without `PS` as
    |will be output as-is.  If `PS` is found in unphased variants, it will be removed.
    |
    |The `--keep-original` option may be used to store the original `PS` value in a new `OPS` field.  The type
    |described in the header will match the original.
  """)
class FixVcfPhaseSet
( @arg(flag='i', doc="Input VCF.") val input: PathToVcf,
  @arg(flag='o', doc="Output VCFs") val output: PathToVcf,
  @arg(flag='s', doc="Samples to mix. See general usage for format and examples.", minElements=0) val samples: Seq[String] = Seq.empty,
  @arg(flag='k', doc="Store the original phase set in the `OPS` field.") val keepOriginal: Boolean = false
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val reader   = VcfSource(input, header=Some(headerForReader))
    val writer   = VcfWriter(output, header=headerForWriter(header=reader.header))
    val progress = ProgressLogger(logger=logger, noun="variants", verb="read", unit=1000000)
    val updater  = new VcfPhaseSetUpdater(header=reader.header, keepOriginal=keepOriginal)

    reader.foreach { variant =>
      progress.record(variant)
      writer.write(variant=updater.update(variant=variant))
    }
    progress.logLast()

    reader.safelyClose
    writer.close()

    logger.info(f"Examined ${updater.resultCounter.total}%,d genotypes.")
    updater.resultCounter.foreach { case (result, count) =>
      logger.info(f"Wrote $count%,d genotypes that were ${result.description}")
    }
  }

  /** Returns a modified [[VcfHeader]] to use when reading the input.  The input VCF header will be modified such
    * that the PS FORMAT field is read as a single fixed value of type [[String]].  This is important, for example, when
    * the header has type [[VcfFieldType.Integer]] but the records have type [[VcfFieldType.String]].
    * */
  def headerForReader: VcfHeader = {
    val reader = VcfSource(input)
    val header = reader.header
    reader.safelyClose()

    val formats = header.formats.filterNot(_.id == "PS")
    val ps      = header.formats
      .find(value => value.id == "PS" && value.count == VcfCount.Fixed(1) && value.kind == VcfFieldType.String)
      .getOrElse {
        VcfFormatHeader(
          id          = "PS",
          count       = VcfCount.Fixed(1),
          kind        = VcfFieldType.String,
          description = "Phasing set (typically the position of the first variant in the set)"
        )
      }
    header.copy(formats=formats :+ ps)
  }

  /** Builds the header for the [[VcfWriter]] from the input reader [[VcfHeader]].  If `keepOriginal` is `true`, then
    * the `OPS` FORMAT header line is added to store the original phase set, which will be a single fixed count of type
    * `String`.  The output phase set will always be a single fixed count of type `Integer`. */
  def headerForWriter(header: VcfHeader):  VcfHeader = {
    // Add the original phase set header FORMAT line if we want to keep the original phase set value
    val original: Option[VcfFormatHeader] = if (!keepOriginal) None else {
      Some(VcfFormatHeader(
        id          = "OPS",
        count       = VcfCount.Fixed(1),
        kind        = VcfFieldType.String,
        description = "Original phasing set (typically the position of the first variant in the set)"
      ))
    }

    // Build the new phase set header FORMAT line.  If the header already contains one that is valid, just keep it.
    val newPhaseSet =  header
      .formats
      .find(value => value.id == "PS" && value.count == VcfCount.Fixed(1) && value.kind == VcfFieldType.Integer)
      .getOrElse {
          VcfFormatHeader(
          id          = "PS",
          count       = VcfCount.Fixed(1),
          kind        = VcfFieldType.Integer,
          description = "Phasing set (typically the position of the first variant in the set)"
        )
      }

    // Put it all together
    header.copy(
      formats = header.formats.filter(_.id != "PS") ++ original.iterator.toSeq :+ newPhaseSet
    )
  }
}

object VcfPhaseSetUpdater {
  /** Base traits for the result of  [[VcfPhaseSetUpdater.updateGenotype()]] */
  sealed trait Result extends EnumEntry {
    def genotype: Genotype
    def description: String
  }
  object Result extends FgBioEnum[Result] {
    case class Valid          (genotype: Genotype) extends Result { val description: String = "was valid" }
    case class Updated        (genotype: Genotype) extends Result { val description: String = "updated" }
    case class MissingPhaseSet(genotype: Genotype) extends Result { val description: String = "phased but missing a phase set" }
    case class MissingPhased  (genotype: Genotype) extends Result { val description: String = "not phased but had a phase set" }
    case class Unphased       (genotype: Genotype) extends Result { val description: String = "not phased" }
    override def values: immutable.IndexedSeq[Result] = findValues
  }
}

class VcfPhaseSetUpdater(header: VcfHeader, keepOriginal: Boolean) extends LazyLogging {
  import VcfPhaseSetUpdater._
  import VcfPhaseSetUpdater.Result._

  private val phaseSetToPositionBySample: Map[String, mutable.HashMap[String, Int]] = {
    header.samples.map { sample => (sample, scala.collection.mutable.HashMap[String, Int]()) }.toMap
  }
  private var lastChrom: String = ""
  private var variants: Long    = 0
  private var genotypes: Long   = 0
  val resultCounter: SimpleCounter[Result] = new SimpleCounter()

  def update(variant: Variant): Variant = {
    // Phase sets cannot span contigs!  So empty our mappings if we get to a new one
    if (variant.chrom != lastChrom) {
      this.phaseSetToPositionBySample.values.foreach(_.clear())
      this.lastChrom = variant.chrom
    }

    variants += 1
    val genotypes: Map[String, Genotype] = variant
      .genotypes
      .map { case (sampleName: String, genotype: Genotype) =>
        // update the genotype
        val result = updateGenotype(variant=variant, genotype=genotype)
        // log the results
        result match {
          case MissingPhased(_)   =>
            logger.warning(f"Genotype had a phase set but was unphased: ${variant.chrom}:${variant.pos}:${variant.id} PS=${genotype("PS")}")
          case MissingPhaseSet(_) =>
            logger.warning(f"Genotype had no phase set but was phased: ${variant.chrom}:${variant.pos}:${variant.id}")
          case _                  => () // do nothing
        }
        this.resultCounter.count(result)
        // return the (potentially) new/updated genotype
        (sampleName, result.genotype)
      }
    variant.copy(genotypes=genotypes)
  }

  private def updateGenotype(variant: Variant, genotype: Genotype): Result = {
    this.genotypes += 1

    (genotype.phased, genotype.get[String]("PS")) match {
      case (true, Some(oldPhaseSet)) => // phased and a phase set
        // Get the new phase set
        val newPhaseSet: Int = this.phaseSetToPositionBySample(genotype.sample).getOrElse(oldPhaseSet, {
          this.phaseSetToPositionBySample(genotype.sample)(oldPhaseSet) = variant.pos
          variant.pos
        })
        // Build the new set of attributes
        val phaseSetAttrs = {
          if (keepOriginal) Map("PS" -> newPhaseSet, "OPS" -> oldPhaseSet)
          else Map("PS" -> newPhaseSet)
        }
        val newGenotype = genotype.copy(attrs = genotype.attrs.filterNot(_._1 == "PS") ++ phaseSetAttrs)
        if (oldPhaseSet == newPhaseSet.toString) Valid(newGenotype) else Updated(newGenotype)
      case (false, Some(ps)) => // warn, not phased but had a phase set!
        MissingPhased(genotype.copy(attrs=genotype.attrs.filterNot(_._1 == "PS"))) // remove the phase set for good measure
      case (true, None) => MissingPhaseSet(genotype) // warn, phased but no phase set
      case (false, None) =>  Unphased(genotype) // ignore, not phased and no phase set

    }
  }
}