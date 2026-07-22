package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.sopt.{arg, clp}
import htsjdk.samtools.util.{Interval, IntervalList}

@clp(group=ClpGroups.Personal, description=
  """
    |Splits intervals within an interval list.  Consumes an input interval list and attempts to break it
    |into _at least_ `--num-intervals` intervals.  The resulting interval list will likely end up with additional
    |intervals.
    |
    |Calculates a max interval size by summing the length of all input intervals and dividing by the desired
    |number of output intervals.  Then each input interval is is broken up into approximately equals size,
    |non-overlapping, sub-intervals of no more than the max size.
  """)
class SplitIntervals(
  @arg(flag='i', doc="Input interval list.") val input: PathToIntervals,
  @arg(flag='o', doc="Output interval list.") val output: PathToIntervals,
  @arg(flag='n', doc="Desired number of output intervals (approx).") val numIntervals: Int)
extends FgBioTool {

  override def execute(): Unit = {
    val intervalList = IntervalList.fromPath(input).uniqued(false)
    val territory = intervalList.getIntervals.map(_.length().toLong).sum
    val maxIntervalSize = Math.ceil(territory / numIntervals.toDouble)

    val outList = new IntervalList(intervalList.getHeader)
    intervalList.getIntervals.foreach { src =>
      if (src.length() <= maxIntervalSize) {
        outList.add(src)
      }
      else {
        val numOutput = Math.ceil(src.length() / maxIntervalSize).toInt
        val eachLength = Math.ceil(src.length() / numOutput.toDouble).toInt

        Range(0, numOutput).map(n => src.getStart + (eachLength * n)).foreach { start =>
          val end = Math.min(start + eachLength - 1, src.getEnd)
          val dest = new Interval(src.getContig, start, end, src.isNegativeStrand, src.getName)
          outList.add(dest)
        }
      }
    }

    val outputTerritory = outList.getIntervals.map(_.length().toLong).sum
    outList.write(output)
    logger.info(f"Input list contained ${intervalList.size()} intervals over ${territory}%,dbp.")
    logger.info(f"Output list contains ${outList.size()} intervals over ${outputTerritory}%,dbp.")
  }
}

