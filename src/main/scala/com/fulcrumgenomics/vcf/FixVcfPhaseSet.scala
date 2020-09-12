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
    |Adds/fixes the phase set (PS) genotype field.
    |
    |The VCF specification allows phased genotypes to be annotated with the `PS` (phase set) `FORMAT` field.  The value
    |should be a non-negative integer, corresponding to the position of the first variant in the phase set.  Some tools
    |will output a non-integer value, as well as describe this field as having non-Integer type in the VCF header.  This
    |tool will update the phase set (`PS`) `FORMAT` field to be VCF spec-compliant.
    |
    |The handling of unphased genotypes with phase-sets is controlled by the `-x` option:
    |- If `-x` is used the genotype will be converted to a phased genotype (e.g. `0/1` => `0|1`)
    |- Otherwise the phase-set (`PS`) will be removed from the genotype
    |
    |The `--keep-original` option may be used to store the original `PS` value in a new `OPS` field.  The type
    |described in the header will match the original.
    |
    |This tool cannot fix phased variants without a phase set, or phased variant sets who have different phase set
    |values.
    |
    |In some cases, VCFs (e.g. from GIAB/NIST or Platinum Genomes) have illegal header lines, for example, a `PEDIGREE`
    |header line without a `ID` key-value field.  The `-z` option can be used to remove those lines.  This option is
    |included in this tool for convenience as those example VCFs in some cases have these illegal header lines, and it
    |is convenient to fix the phase set in addition to removing those illegal header lines.
  """)
class FixVcfPhaseSet
( @arg(flag='i', doc="Input VCF.") val input: PathToVcf,
  @arg(flag='o', doc="Output VCF.") val output: PathToVcf,
  @arg(flag='k', doc="Store the original phase set in the `OPS` field.") val keepOriginal: Boolean = false,
  @arg(flag='x', doc="Set unphased genotypes with a PS FORMAT value to be phased.") val phaseGenotypesWithPhaseSet: Boolean = false,
  @arg(flag='z', doc="Remove header lines that do not contain an ID key-value, which is required in VCF.") val removeNoIdHeaderLines: Boolean = false
) extends FgBioTool with LazyLogging {

  import FixVcfPhaseSet._

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    // Get the original VCF header
    val inHeader = {
      val source = VcfSource(input)
      val header = source.header
      source.safelyClose()
      header
    }
    // Developer note: we modify the input header so we are able to read in the PS format field as a string.
    val reader   = VcfSource(input, headerTransformer=headerForReader)
    val writer   = VcfWriter(output, header=headerForWriter(header=inHeader))
    val progress = ProgressLogger(logger=logger, noun="variants", verb="read", unit=100000)
    val updater  = new VcfPhaseSetUpdater(header=reader.header, keepOriginal=keepOriginal, phaseGenotypesWithPhaseSet=phaseGenotypesWithPhaseSet)

    reader.foreach { variant =>
      writer.write(variant=updater.update(variant=variant))
      progress.record(variant)
    }
    progress.logLast()

    reader.safelyClose()
    writer.close()

    logger.info(f"Examined ${updater.descriptionCounter.total}%,d genotypes.")
    updater.descriptionCounter.foreach { case (description, count) =>
      logger.info(f"Wrote $count%,d genotypes that $description")
    }
  }

  /** Returns a modified [[VcfHeader]] to use when reading the input.  The input VCF header will be modified such
    * that the PS FORMAT field is read as a single fixed value of type [[String]].  This is important, for example, when
    * the header has type [[VcfFieldType.Integer]] but the records have type [[VcfFieldType.String]].
    * */
  def headerForReader(header: VcfHeader): VcfHeader = {
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
        description = "Original phasing set prior to being fixed with fgbio FixVcfPhaseSet"
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

    // In some cases, input VCFs may have header lines (ex. PEDIGREE) with no ID, which is required
    val others = if (removeNoIdHeaderLines) header.others.filterNot(_.id == null) else header.others

    // Put it all together
    header.copy(
      formats = header.formats.filter(_.id != "PS") ++ original.toSeq :+ newPhaseSet,
      others  = others
    )
  }
}

object FixVcfPhaseSet {
  private[vcf] object VcfPhaseSetUpdater {
    /** Base traits for the result of  [[VcfPhaseSetUpdater.updateGenotype()]] */
    sealed trait Result extends EnumEntry {
      def genotype: Genotype
      def description: String
    }
    object Result extends FgBioEnum[Result] {
      case class Valid                      (genotype: Genotype) extends Result { val description: String = "were valid" }
      case class UpdatedPhaseSetValue       (genotype: Genotype) extends Result { val description: String = "had phase set value updated only" }
      case class UpdatedPhaseAndValue       (genotype: Genotype) extends Result { val description: String = "were set to phased and had phase set value updated" }
      case class UpdatedPhaseOnly           (genotype: Genotype) extends Result { val description: String = "were set to phased" }
      case class PhasedMissingPhaseSetValue (genotype: Genotype) extends Result { val description: String = "were phased but missing a phase set value" }
      case class NotPhasedWithPhaseSetValue (genotype: Genotype) extends Result { val description: String = "were not phased but had a phase set value" }
      case class NotPhasedNoPhaseSetValue   (genotype: Genotype) extends Result { val description: String = "were not phased and no phase set value" }
      override def values: immutable.IndexedSeq[Result] = findValues
    }
  }

  /** Updates the phase set for a given variant's genotype(s).
    *
    * The phase set values must have type [[String]].
    *
    * @param header the VCF header
    * @param keepOriginal true to store the original phase set value in the OPS format field, false otherwise
    * @param phaseGenotypesWithPhaseSet set unphased genotypes with a PS FORMAT value to be phased
    */
  private[vcf] class VcfPhaseSetUpdater(header: VcfHeader, keepOriginal: Boolean, phaseGenotypesWithPhaseSet: Boolean) extends LazyLogging {
    import VcfPhaseSetUpdater._
    import VcfPhaseSetUpdater.Result._

    private val phaseSetToPositionBySample: Map[String, mutable.HashMap[String, Int]] = {
      header.samples.map { sample => (sample, scala.collection.mutable.HashMap[String, Int]()) }.toMap
    }
    private var lastChrom: String = ""
    private var lastPos: Int      = 0
    private var variants: Long    = 0
    private var genotypes: Long   = 0
    val descriptionCounter: SimpleCounter[String] = new SimpleCounter()

    /** Updates the phase set of the variant's genotype(s) */
    def update(variant: Variant): Variant = {
      // Phase sets cannot span contigs!  So empty our mappings if we get to a new one
      if (variant.chrom != lastChrom) {
        this.phaseSetToPositionBySample.values.foreach(_.clear())
        this.lastChrom = variant.chrom
        this.lastPos   = variant.pos
      }
      assert(this.lastPos <= variant.pos, "Variants are out of order!")

      variants += 1
      val genotypes: Map[String, Genotype] = variant.genotypes
        .map { case (sampleName: String, genotype: Genotype) =>
          // update the genotype
          val result = updateGenotype(variant=variant, genotype=genotype)
          // log the results
          result match {
            case NotPhasedWithPhaseSetValue(_)   =>
              logger.debug(
                "Genotype had a phase set but was unphased:" +
                  f"${variant.chrom}:${variant.pos}:${variant.id.getOrElse(".")} sample=$sampleName PS=${genotype("PS")}" +
                  "; consider r-running using `-x/--phase-genotypes-with-phase-set`"
              )
            case PhasedMissingPhaseSetValue(_) =>
              logger.debug(
                f"Genotype had no phase set but was phased: ${variant.chrom}:${variant.pos}:${variant.id.getOrElse(".")}"
              )
            case _                  => () // do nothing
          }
          this.descriptionCounter.count(result.description)
          // return the (potentially) new/updated genotype
          (sampleName, result.genotype)
        }
      variant.copy(genotypes=genotypes)
    }

    /** Updates the phase set of a given genotype.  Returns the [[Result]] of updating the genotyping, describing
      * what updating was done, if any. */
    private[vcf] def updateGenotype(variant: Variant, genotype: Genotype): Result = {
      this.genotypes += 1

      (genotype.phased, genotype.get[String]("PS")) match {
        case (true, None)  => PhasedMissingPhaseSetValue(genotype) // phased but no phase set, what can we do?
        case (false, None) => NotPhasedNoPhaseSetValue(genotype)   // not phased and no phase set, everything is fine
        case (false, Some(_)) if !phaseGenotypesWithPhaseSet =>    // unphased and a phase set, but we do not want to convert it to phased
          NotPhasedWithPhaseSetValue(genotype.copy(attrs=genotype.attrs.filterNot(_._1 == "PS"))) // remove the phase set for good measure
        case (isPhased, Some(oldValue)) => // may or not be phased, but has a phase set
          // Get the new phase set
          val phaseSetToValue = this.phaseSetToPositionBySample(genotype.sample)
          val newValue        = phaseSetToValue.getOrElseUpdate(oldValue, variant.pos)
          // Build the new set of attributes
          val phaseSetAttrs = {
            if (keepOriginal) Map("PS" -> newValue, "OPS" -> oldValue)
            else Map("PS" -> newValue)
          }
          // update the genotype
          val newGenotype = genotype.copy(attrs = genotype.attrs.filterNot(_._1 == "PS") ++ phaseSetAttrs, phased=true)
          // return based on if we phased the variant and if we added/updated the phase set value
          (isPhased, oldValue == newValue.toString) match {
            case (true, true)   => Valid(newGenotype)
            case (true, false)  => UpdatedPhaseSetValue(newGenotype)
            case (false, true)  => UpdatedPhaseOnly(newGenotype)
            case (false, false) => UpdatedPhaseAndValue(newGenotype)
          }
      }
    }
  }
}
