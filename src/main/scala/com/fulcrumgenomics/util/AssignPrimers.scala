package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToBam, SafelyClosable}
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}


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
    |Each read is assigned independently of its mate (for paired end reads). The primer for a read is assumed to be
    |located at the start of the read in 5' sequencing order.  Therefore, a positive strand
    |read will use its aligned start position to match against the amplicon's left-most coordinate, while a negative
    |strand read will use its aligned end position to match against the amplicon's right-most coordinate.
    |
    |For paired end reads, the assignment for mate will also be stored in the current read, using the same procedure as
    |above but using the mate's coordinates.  This requires the input BAM have the mate-cigar ("MC") SAM tag.  Read
    |pairs must have both ends mapped in forward/reverse configuration to have an assignment.  Furthermore, the amplicon
    |assignment may be different for a read and its mate.  This may occur, for example, if tiling nearby amplicons and
    |a large deletion occurs over a given primer and therefore "skipping" an amplicon.  This may also occur if there are
    |translocations across amplicons.
    |
    |The output will have the following tags added:
    |- ap: the assigned primer coordinates (ex. `chr1:1010873-1010894`)
    |- am: the mate's assigned primer coordinates (ex. `chr1:1011118-1011137`)
    |- ip: the assigned amplicon id
    |- im: the mate's assigned amplicon id (or `=` if the same as the assigned amplicon)
    |
    |The read sequence of the primer is not checked against the expected reference sequence at the primer's genomic
    |coordinates.
    |
    |In some cases, large deletions within one end of a read pair may cause a primary and supplementary alignments to be
    |produced by the aligner, with the supplementary alignment containing the primer end of the read (5' sequencing order).
    |In this case, the primer may not be assigned for this end of the read pair.  Therefore, it is recommended to prefer
    |or choose the primary alignment that has the closest aligned read base to the 5' end of the read in sequencing order.
    |For example, from `bwa` version `0.7.16` onwards, the `-5` option may be used.  Consider also using the `-q` option 
    |for `bwa` `0.7.16` as well, which is standard in `0.7.17` onwards when the `-5` option is used.
    |
    |The `--annotate-all` option may be used to annotate all alignments for a given read end (eg. R1) with
    |the same assignment.  If the assignment differs across alignments for the same read end, no assignment is given.
    |Furthermore, if the input BAM is neither `queryname` sorted nor `query` grouped, it will be sorted into queryname
    |order to assign all alignments cross a template simultaneously.  The output is written in coordinate order.
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
 @arg(doc="Annotate all R1 (or R2) with same value.") val annotateAll: Boolean = false
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
      detector             = Amplicon.overlapDetector(amplicons=amplicons.iterator),
      slop                 = slop,
      unclippedCoordinates = unclippedCoordinates
    )

    val labeller = new AmpliconLabeller(
      amplicons                 = amplicons.iterator,
      primerCoordinatesTag      = primerCoordinatesTag,
      matePrimerCoordinatesTag  = matePrimerCoordinatesTag,
      ampliconIdentifierTag     = ampliconIdentifierTag,
      mateAmpliconIdentifierTag = mateAmpliconIdentifierTag
    )

    if (annotateAll) {
      // do not sort if the input order is query name or query grouped or unsorted
      val sorter = SamOrder(reader.header).flatMap { order =>
          if (order.sortOrder == SortOrder.queryname) None
          else if (order.groupOrder == GroupOrder.query ) None
          else if (order == SamOrder.Unsorted) None
          else Some(Bams.sorter(order, reader.header))
      }
      val writer = SamWriter(output, reader.header)
      // All alignments will be assigned a primer independently, and then grouped by read end. An primer assignment will
      // be made on a given read end only if there is one and only one unique primer assignment.
      Bams.templateIterator(reader).foreach { template =>
        // split into R1s and R2s
        val (r1s, r2s) = template.allReads.toSeq.partition(r => r.unpaired || r.firstOfPair)
        // detect the primers for each of the R1s and R2s
        val r1Amplicon  = r1s.flatMap(detector.findPrimer).distinct match {
          case Seq(x) => Some(x)
          case _      => None
        }
        val r2Amplicon  = r2s.flatMap(detector.findPrimer).distinct match {
          case Seq(x) => Some(x)
          case _      => None
        }
        // label the reads
        r1s.foreach { rec => labeller.label(rec=rec, recAmplicon=r1Amplicon, mateAmplicon=r2Amplicon) }
        r2s.foreach { rec => labeller.label(rec=rec, recAmplicon=r2Amplicon, mateAmplicon=r1Amplicon) }
        // write them out
        template.allReads.foreach { r =>
          sorter match {
            case Some(_sorter) => _sorter += r
            case None          => writer += r
          }
          progress.record(r)
        }
      }
      progress.logLast()
      // Write them out after sorting, if necessary
      sorter.foreach { _sorter => writer ++= _sorter }
      writer.close()
    }
    else {
      val writer = SamWriter(output, reader.header)

      reader.foreach { rec =>
        val recAmplicon = detector.findPrimer(rec=rec)
        labeller.label(
          rec          = rec,
          recAmplicon  = recAmplicon,
          mateAmplicon = if (rec.unpaired) None else detector.findMatePrimer(rec=rec)
        )
        writer += rec
        progress.record(rec)
      }
      progress.logLast()
      writer.close()
    }

    reader.safelyClose()

    // Sum up the values across all amplicons
    val totalMetric = new AssignPrimersMetric(identifier=AssignPrimersMetric.AllAmpliconsIdentifier)
    labeller.metrics.values.foreach { metric =>
      totalMetric.r1s   += metric.r1s
      totalMetric.r2s   += metric.r2s
      totalMetric.left  += metric.left
      totalMetric.right += metric.right
      totalMetric.pairs += metric.pairs
    }

    // Write the metrics
    Metric.write[AssignPrimersMetric](
      path    = metrics,
      metrics = (amplicons.iterator.map(labeller.metrics.apply) ++ Iterator(totalMetric)).map(_.finalize(total=progress.getCount))
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
  val PrimerCoordinateTag      : String = "rp"
  /** The SAM tag to use for the current record's mate's assigned primer genomic coordinates */
  val MatePrimerCoordinateTag  : String = "mp"
  /** The SAM tag to use for the current record's assigned primer identifier */
  val AmpliconIdentifierTag    : String = "ra"
  /** The SAM tag to use for the current record's mate's assigned primer identifier */
  val MateAmpliconIdentifierTag: String = "ma"

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
  * @param amplicons the amplicons used to collect metrics, may be empty
  * @param primerCoordinatesTag the SAM tag to store the primer genomic coordinates for the read
  * @param matePrimerCoordinatesTag the SAM tag to store the primer genomic coordinates for the read's mate
  * @param ampliconIdentifierTag the SAM tag to store the amplicon identifier for the read
  * @param mateAmpliconIdentifierTag the SAM tag to store  the amplicon identifier for the read's mate
  */
class AmpliconLabeller(val amplicons: Iterator[Amplicon] = Iterator.empty,
                       val primerCoordinatesTag: String = AssignPrimers.PrimerCoordinateTag,
                       val matePrimerCoordinatesTag: String = AssignPrimers.MatePrimerCoordinateTag,
                       val ampliconIdentifierTag: String = AssignPrimers.AmpliconIdentifierTag,
                       val mateAmpliconIdentifierTag: String = AssignPrimers.MateAmpliconIdentifierTag,
                      ) {

  val metrics: Map[Amplicon, AssignPrimersMetric] = amplicons.map { amplicon =>
    amplicon -> new AssignPrimersMetric(identifier=amplicon.identifier)
  }.toMap

  /** Labels (adds SAM tags) to a record based on the assigned amplicon(s) for itself and potentially its mate.
    *
    * @param rec the record to label
    * @param recAmplicon the assigned amplicon for the record
    * @param mateAmplicon the assigned amplicon fo the record's mate (if paired)
    */
  def label(rec: SamRecord,
            recAmplicon: Option[Amplicon] = None,
            mateAmplicon: Option[Amplicon] = None): Unit = {
    // Process the primer for the current read
    recAmplicon.foreach { amp =>
      // update metrics
      if (rec.primary && !rec.supplementary) {
        metrics.get(amp).foreach { metric =>
          if (rec.positiveStrand) metric.left += 1 else metric.right += 1
          if (rec.unpaired || rec.firstOfPair) metric.r1s += 1 else metric.r2s += 1
        }
      }
      // assign ap/ip
      rec(primerCoordinatesTag)  = if (rec.positiveStrand) amp.leftPrimerLocation else amp.rightPrimerLocation
      rec(ampliconIdentifierTag) = amp.identifier
    }

    // Find the primer for its mate
    mateAmplicon.foreach { amp =>
      require(rec.paired)
      val isPrimerPair = rec.isFrPair && recAmplicon.contains(amp)
      // update metrics
      metrics.get(amp).foreach { metric =>
        if (isPrimerPair && rec.firstOfPair) metric.pairs += 1
      }
      // assign am/im
      rec(matePrimerCoordinatesTag)   = if (rec.matePositiveStrand) amp.leftPrimerLocation else amp.rightPrimerLocation
      rec(mateAmpliconIdentifierTag) = if (isPrimerPair) "=" else amp.identifier
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

/** Metrics produced by `AssignPrimers` that detail how many reads were assigned to a given primer and/or amplicon.
  *
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
  var frac_left: Double = 0,
  var frac_right: Double = 0,
  var frac_r1s: Double = 0,
  var frac_r2s: Double = 0,
  var frac_pairs: Double = 0
) extends Metric {
  /** Update the factional metrics given the total number of reads. */
  def finalize(total: Long): this.type = {
    if (total == 0) {
      this.frac_left  = 0
      this.frac_right = 0
      this.frac_r1s   = 0
      this.frac_r2s   = 0
      this.frac_pairs = 0
    }
    else {
      this.frac_left  = this.left / total.toDouble
      this.frac_right = this.right / total.toDouble
      this.frac_r1s   = this.r1s / total.toDouble
      this.frac_r2s   = this.r2s / total.toDouble
      this.frac_pairs = if (total <= 1) 0 else this.pairs / (total / 2.0)
    }
    this
  }
}

object AssignPrimersMetric {
  /** The name to use for the [[AssignPrimersMetric]] calculated over all amplicons */
  val AllAmpliconsIdentifier: String = "AllAmplicons"
}
