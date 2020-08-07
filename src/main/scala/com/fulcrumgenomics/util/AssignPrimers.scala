package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}


@clp(group = ClpGroups.SamOrBam, description=
  """
    |Assigns reads to primers post-alignment. Takes in a BAM file of aligned reads and a tab-delimited file with five columns
    |(`chrom`, `left_start`, `left_end`, `right_start`, and `right_end`) which provide the 1-based inclusive start and
    |end positions of the primers for each amplicon.  The primer file must include headers, e.g:
    |
    |```
    |chrom  left_start  left_end  right_start right_end
    |chr1   1010873     1010894   1011118     1011137
    |```
    |
    |Optionally, a sixth column column `id` may be given with a unique name for the amplicon.  If not given, the
    |coordinates of the amplicon's primers will be used:
    |  `<chrom>:<left_start>-<left_end>,<chrom>:<right_start>:<right_end>`
    |
    |The output will have the following tags added:
    |- ap: the assigned primer coordinates (ex. `chr1:1010873-1010894`)
    |- am: the mate's assigned primer coordinates (ex. `chr1:1011118-1011137`)
    |- ip: the assigned amplicon id
    |- im: the mate's assigned amplicon id (or `=` if the same as the assigned amplicon)
    |
    |The read sequence of the primer is not checked against the expected reference sequence at the primer's genomic
    |coordinates.
  """)
class AssignPrimers
(@arg(flag='i', doc="Input BAM file.")  val input: PathToBam,
 @arg(flag='o', doc="Output BAM file.") val output: PathToBam,
 @arg(flag='m', doc="Output metrics file.") val metrics: FilePath,
 @arg(flag='p', doc="File with primer locations.") val primers: FilePath,
 @arg(flag='S', doc="Match to primer locations +/- this many bases.") val slop: Int = 5,
 @arg(flag='U', doc="True to based on the unclipped coordinates (adjust based on hard/soft clipping), otherwise the aligned bases") val unclippedCoordinates: Boolean = true,
 @arg(doc="The SAM tag for the assigned primer coordinate.") val primerCoordinatesTag: String = AssignPrimers.PrimerCoordinateTag,
 @arg(doc="The SAM tag for the mate's assigned primer coordinate.") val matePrimerCoordinatesTag: String = AssignPrimers.MatePrimerCoordinateTag,
 @arg(doc="The SAM tag for the assigned amplicon identifier.") val ampliconIdentifierTag: String = AssignPrimers.AmpliconIdentifierTag,
 @arg(doc="The SAM tag for the mate's assigned amplicon identifier.") val mateAmpliconIdentifierTag: String = AssignPrimers.MateAmpliconIdentifierTag,
) extends FgBioTool with LazyLogging {

  Io.assertReadable(input)
  Io.assertReadable(primers)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val reader    = SamSource(input)
    val writer    = SamWriter(output, reader.header)
    val progress  = ProgressLogger(logger=logger, unit=250000)
    val amplicons = Metric.read[Amplicon](path=primers)
    val detector  = new AmpliconDetector(
      detector             = Amplicon.detector(amplicons=amplicons.iterator),
      slop                 = slop,
      unclippedCoordinates = unclippedCoordinates
    )

    val metricsMap: Map[Amplicon, AssignPrimersMetric] = amplicons.map { amplicon =>
      amplicon -> new AssignPrimersMetric(identifier=amplicon.identifier)
    }.toMap

    val labeller = new AmpliconLabeller(
      metrics                   = metricsMap,
      primerCoordinatesTag      = primerCoordinatesTag,
      matePrimerCoordinatesTag  = matePrimerCoordinatesTag,
      ampliconIdentifierTag     = ampliconIdentifierTag,
      mateAmpliconIdentifierTag = mateAmpliconIdentifierTag
    )

    reader.foreach { rec =>
      val recAmplicon = detector.find(rec=rec)
      labeller.label(
        rec = rec,
        recAmplicon = recAmplicon,
        mateAmplicon = if (rec.unpaired) None else detector.findMate(rec=rec)
      )
      writer.write(rec)
      progress.record(rec)
    }
    progress.logLast()

    reader.safelyClose()
    writer.close()

    // Sum up the values across all amplicons
    val totalMetric = new AssignPrimersMetric(identifier=AssignPrimersMetric.AllAmpliconsIdentifier)
    metricsMap.values.foreach { metric =>
      totalMetric.r1s   += metric.r1s
      totalMetric.r2s   += metric.r2s
      totalMetric.left  += metric.left
      totalMetric.right += metric.right
      totalMetric.pairs += metric.pairs
    }

    // Write the metrics
    Metric.write[AssignPrimersMetric](
      path    = metrics,
      metrics = (amplicons.iterator.map(metricsMap.apply) ++ Iterator(totalMetric)).map(_.finalize(total=progress.getCount))
    )

    // Log some info
    def log(numerator: Long, noun: String): Unit = {
      val pct: Double = if (progress.getCount == 0) 0 else numerator * 100.0 / progress.getCount
      logger.info(f"Assigned $numerator%,d out of ${progress.getCount}%,d ($pct%.2f%%) reads $noun.")
    }
    log(numerator=totalMetric.left, "to left primers")
    log(numerator=totalMetric.right, "to right primers")
    log(numerator=totalMetric.pairs * 2, "as primer pairs")
  }
}

object AssignPrimers {
  /** The SAM tag to use for the current record's assigned primer genomic coordinates */
  val PrimerCoordinateTag      : String = "ap"
  /** The SAM tag to use for the current record's mate's assigned primer genomic coordinates */
  val MatePrimerCoordinateTag  : String = "am"
  /** The SAM tag to use for the current record's assigned primer identifier */
  val AmpliconIdentifierTag    : String = "ip"
  /** The SAM tag to use for the current record's mate's assigned primer identifier */
  val MateAmpliconIdentifierTag: String = "im"

  def tags: Seq[String] = Seq(PrimerCoordinateTag, MatePrimerCoordinateTag, AmpliconIdentifierTag, MateAmpliconIdentifierTag)

}

/** Labels (add SAM tags) to records based on assigned amplicons.
  *
  * The output will have the following tags added:
  * - ap: the assigned primer coordinates (ex. `chr1:1010873-1010894`)
  * - am: the mate's assigned primer coordinates (ex. `chr1:1011118-1011137`)
  * - ip: the assigned amplicon id
  * - im: the mate's assigned amplicon id (or `=` if the same as the assigned amplicon)
  *
  * If not primer/amplicon was found, then no tag will be written.
  *
  * @param metrics a map of amplicon to [[AssignPrimersMetric]] to collect metrics, may be empty
  * @param primerCoordinatesTag the SAM tag to store the primer genomic coordinates for the read
  * @param matePrimerCoordinatesTag the SAM tag to store the primer genomic coordinates for the read's mate
  * @param ampliconIdentifierTag the SAM tag to store the amplicon identifier for the read
  * @param mateAmpliconIdentifierTag the SAM tag to store  the amplicon identifier for the read's mate
  */
class AmpliconLabeller(val metrics: Map[Amplicon, AssignPrimersMetric] = Map.empty,
                       val primerCoordinatesTag: String = AssignPrimers.PrimerCoordinateTag,
                       val matePrimerCoordinatesTag: String = AssignPrimers.MatePrimerCoordinateTag,
                       val ampliconIdentifierTag: String = AssignPrimers.AmpliconIdentifierTag,
                       val mateAmpliconIdentifierTag: String = AssignPrimers.MateAmpliconIdentifierTag,
                      ) {
  /** Labels (adds SAM tags) to a record based on the assigned amplicon(s) for itself and potentially its mate.
    *
    * @param rec the record to label
    * @param recAmplicon the assigned amplicon for the record
    * @param mateAmplicon the assigned amplicon fo the record's mate (if paired)
    */
  def label(rec: SamRecord,
            recAmplicon: Option[Amplicon] = None,
            mateAmplicon: Option[Amplicon] = None): Unit = {
    // Find the primer for the current read
    recAmplicon.foreach { amp =>
      // update metrics
      metrics.get(amp).foreach { metric =>
        if (rec.positiveStrand) metric.left += 1 else metric.right += 1
        if (rec.unpaired || rec.firstOfPair) metric.r1s += 1 else metric.r2s += 1
      }
      // assign ap/ip
      rec(primerCoordinatesTag)   = if (rec.positiveStrand) amp.leftPrimerString else amp.rightPrimerString
      rec(ampliconIdentifierTag) = amp.identifier
    }

    if (rec.paired) {
      // Find the primer for its mate
      mateAmplicon.foreach { amp =>
        val isPrimerPair = rec.isFrPair && recAmplicon.contains(amp)
        // update metrics
        metrics.get(amp).foreach { metric =>
          if (isPrimerPair && rec.firstOfPair) metric.pairs += 1
        }
        // assign am/im
        rec(matePrimerCoordinatesTag)   = if (rec.matePositiveStrand) amp.leftPrimerString else amp.rightPrimerString
        rec(mateAmpliconIdentifierTag) = if (isPrimerPair) "=" else amp.identifier
      }
    }
  }

  /** Labels (adds SAM tags) to a read pair based on the assigned amplicon(s) for itself and potentially its mate.
    *
    * @param r1 the first of pair
    * @param r2 the second of pair
    * @param recAmplicon the assigned amplicon for the record
    * @param mateAmplicon the assigned amplicon fo the record's mate (if paired)
    */
  def label(r1: SamRecord,
            r2: SamRecord,
            recAmplicon: Option[Amplicon],
            mateAmplicon: Option[Amplicon]): Unit = {
    label(rec=r1, recAmplicon=recAmplicon, mateAmplicon=mateAmplicon)
    label(rec=r2, recAmplicon=mateAmplicon, mateAmplicon=recAmplicon)
  }
}

/**
  * @param identifier the amplicon identifier this metric collects over
  * @param left the number of reads assigned to the left primer
  * @param right the number of reads assigned to the right primer
  * @param r1s the number of R1 reads assigned to this amplicon
  * @param r2s the number of R2 reads assigned to this amplicon
  * @param pairs the number of read pairs where R1 and R2 are both assigned to the this amplicon and are in FR orientation
  * @param frac_left the fraction of reads assigned to the left primer
  * @param frac_right the fraction of reads assigned to the right primer
  * @param frac_r1s the fraction of R1s reads assigned to this amplicon
  * @param frac_r2s the fraction of R2s reads assigned to this amplicon
  * @param frac_pairs the fraction of read pairs where R1 and R2 are both assigned to the this amplicon and are in FR orientation
  */
case class AssignPrimersMetric
( identifier: String,
  var left: Long = 0L,
  var right: Long = 0L,
  var r1s: Long = 0L,
  var r2s: Long = 0L,
  var pairs: Long = 0L,
  var frac_left: Double = 0L,
  var frac_right: Double = 0L,
  var frac_r1s: Double = 0L,
  var frac_r2s: Double = 0L,
  var frac_pairs: Double = 0L
) extends Metric {
  /** Update the factional metrics given the total number of reads. */
  def finalize(total: Long): this.type = if (total == 0) this else {
    this.frac_left  = this.left / total
    this.frac_right = this.right / total
    this.frac_r1s   = this.r1s / total
    this.frac_r2s   = this.r2s / total
    this.frac_pairs = if (total <= 1) 0 else this.pairs / (total / 2)
    this
  }
}

object AssignPrimersMetric {
  /** The name to use for the [[AssignPrimersMetric]] calculated over all amplicons */
  val AllAmpliconsIdentifier: String = "AllAmplicons"
}

