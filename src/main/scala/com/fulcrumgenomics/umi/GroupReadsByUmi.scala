/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package com.fulcrumgenomics.umi

import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.Bams
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamSource, SamWriter}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.{LazyLogging, NumericCounter, SimpleCounter}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.umi.GroupReadsByUmi._
import com.fulcrumgenomics.util.Metric.{Count, Proportion}
import com.fulcrumgenomics.util.Sequences.countMismatches
import com.fulcrumgenomics.util._
import enumeratum.EnumEntry
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.util.SequenceUtil

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable
import scala.collection.mutable.ListBuffer


object GroupReadsByUmi {
  /** A type representing an actual UMI that is a string of DNA sequence. */
  type Umi = String
  /** A type to represent a unique ID assigned to a molecule. */
  type MoleculeId = String

  /** The suffix appended to molecular identifier reads from the "top strand" of the duplex molecule. These reads have
    * 5' coordinate for read 1 less than or equal to the 5' coordinate for read 2. */
  val TopStrandDuplex = "/A"

  /** The suffix appended to molecular identifier reads from the "bottom strand" of the duplex molecule. These reads have
    * 5' coordinate for read 1 greater than the 5' coordinate for read 2. */
  val BottomStrandDuplex = "/B"

  private val ReadInfoTempAttributeName = "__GRBU_ReadInfo"

  /** A case class to represent all the information we need to order reads for duplicate marking / grouping. */
  case class ReadInfo(refIndex: Int, start1: Int, start2: Int, strand1: Boolean, strand2: Boolean, library: String)

  object ReadInfo {
    /** Looks in all the places the library name can be hiding. Returns the library name
      * if one is found, otherwise returns "unknown".
      */
    private def library(rec: SamRecord): String = {
      val rg = rec.readGroup
      if (rg != null && rg.getLibrary != null) rg.getLibrary else "unknown"
    }

    /** Creates/retrieves a ReadEnds object from a SamRecord and stores it in a temporary attribute for later user. */
    def apply(rec: SamRecord) : ReadInfo = {
      val tmp = rec.transientAttrs[ReadInfo](GroupReadsByUmi.ReadInfoTempAttributeName)
      if (tmp != null) {
        tmp
      }
      else {
        if (rec.refIndex != rec.mateRefIndex) throw new IllegalArgumentException("Mate on different chrom.")
        val chrom   = rec.refIndex
        val recNeg  = rec.negativeStrand
        val recPos  = if (recNeg) rec.unclippedEnd else rec.unclippedStart
        val mateNeg = rec.mateNegativeStrand
        val matePos = if (mateNeg) SAMUtils.getMateUnclippedEnd(rec.asSam) else SAMUtils.getMateUnclippedStart(rec.asSam)
        val lib = library(rec)

        val result = if (recPos < matePos || (recPos == matePos && !recNeg)) {
          new ReadInfo(chrom, recPos, matePos, recNeg, mateNeg, lib)
        }
        else {
          new ReadInfo(chrom, matePos, recPos, mateNeg, recNeg, lib)
        }

        rec.transientAttrs(GroupReadsByUmi.ReadInfoTempAttributeName, result)
        result
      }
    }
  }

  /** Trait that can be implemented to provide a UMI assignment strategy. */
  private[umi] sealed trait UmiAssigner {
    private val counter = new AtomicLong()

    /** Take in a sequence of UMIs and assign each UMI to a unique UMI group ID. */
    def assign(rawUmis: Seq[Umi]) : Map[Umi, MoleculeId]

    /** Default implementation of a method to retrieve the next ID based on a counter. */
    protected def nextId: MoleculeId = this.counter.getAndIncrement().toString

    /** Takes in a map of UMI to "sentinel" UMI, and outputs a map of UMI -> Molecule ID. */
    protected def assignIds(assignments: Map[Umi,Umi]): Map[Umi,MoleculeId] = {
      val idMap = assignments.values.toSet.map((umi:Umi) => umi -> nextId).toMap
      assignments.map { case (k,v) => (k, idMap(v)) }
    }
  }

  /**
    * Assigns UMIs based solely on sequence identity.
    */
  private[umi] class IdentityUmiAssigner extends UmiAssigner {
    override def assign(rawUmis: Seq[Umi]): Map[Umi, MoleculeId] = assignIds(rawUmis.map(umi => umi -> umi.toUpperCase).toMap)
  }

  /**
    * Assigns UMIs based solely on sequence identity.
    */
  private[umi] class SimpleErrorUmiAssigner(val maxMismatches: Int) extends UmiAssigner {
    override def assign(rawUmis: Seq[Umi]): Map[Umi, MoleculeId] = {
      type UmiSet = mutable.SortedSet[Umi]
      val umiSets = ListBuffer[UmiSet]()
      for (umi <- rawUmis) {
        val matchedSets = mutable.Set[UmiSet]()
        for (set <- umiSets) {
          if (set.exists(other => countMismatches(other, umi) <= maxMismatches)) {
            set.add(umi)
            matchedSets.add(set)
          }
        }

        matchedSets.size match {
          case 0 => umiSets += mutable.SortedSet(umi)
          case 1 => matchedSets.head.add(umi)
          case _ =>
            // Merge the sets and remove all but one from umiSets
            val pick = matchedSets.head
            matchedSets.tail.foreach(set => {
              umiSets -= set
              pick ++= set
            })
        }
      }

      assignIds(umiSets.flatMap(set => set.map(umi => umi -> set.head)).toMap)
    }
  }


  /**
    * Class that implements the directed adjacency graph method from umi_tools.
    * See: https://github.com/CGATOxford/UMI-tools
    */
  private[umi] class AdjacencyUmiAssigner(val maxMismatches: Int) extends UmiAssigner {
    /** Represents a node in the adjacency graph; equality is just by UMI sequence. */
    class Node(val umi: Umi, val count: Long, val children: mutable.Buffer[Node] = mutable.Buffer()) {
      /** Gets the full set of descendants from this node. */
      def descendants: List[Node] = {
        val buffer = ListBuffer[Node]()
        addDescendants(buffer)
        buffer.toList
      }

      private def addDescendants(buffer: mutable.Buffer[Node]): Unit = children.foreach(child => {
        buffer += child
        child.addDescendants(buffer)
      })

      override def equals(other: scala.Any): Boolean = other.isInstanceOf[Node] && this.umi == other.asInstanceOf[Node].umi

      override def toString: String = s"$umi [$count]"
    }

    /** Returns the count of each raw UMI that was observed. */
    protected def count(umis: Seq[Umi]): Iterator[(Umi,Long)] = SimpleCounter(umis).iterator

    /** Returns whether or not a pair of UMIs match closely enough to be considered adjacent in the graph. */
    protected def matches(lhs: Umi, rhs: Umi): Boolean = {
      require(lhs.length == rhs.length, s"UMIs of different length detected: $lhs vs. $rhs")
      var idx = 0
      var mismatches = 0
      val len = lhs.length
      while (idx < len && mismatches <= this.maxMismatches) {
        if (lhs(idx) != rhs(idx)) mismatches += 1
        idx += 1
      }

      mismatches <= maxMismatches
    }

    /** Assigns IDs to each UMI based on the root to which is it mapped. */
    protected def assignIdsToNodes(roots: Seq[Node]): Map[Umi, MoleculeId] = {
      val mappings = mutable.Buffer[(Umi,MoleculeId)]()
      roots.foreach(root => {
        val id = nextId
        mappings.append((root.umi, id))
        root.descendants.foreach(child => mappings.append((child.umi, id)))
      })

      mappings.toMap
    }

    override def assign(rawUmis: Seq[Umi]): Map[Umi, MoleculeId] = {
      // A list of all the root UMIs/Nodes that we find
      val roots = mutable.Buffer[Node]()

      // Make a list of counts of all UMIs in order from most to least abundant; we'll consume from this buffer
      var remaining = count(rawUmis).map{ case(umi,count) => new Node(umi, count) }.toBuffer.sortBy((n:Node) => -n.count)

      // Now build one or more graphs starting with the most abundant remaining umi
      while (remaining.nonEmpty) {
        val nextRoot = remaining.remove(0)
        roots += nextRoot
        val working = mutable.Buffer[Node](nextRoot)

        while (working.nonEmpty) {
          val root = working.remove(0)
          val (hits, misses) = remaining.partition(other => root.count >= 2 * other.count - 1 && matches(root.umi, other.umi))
          root.children ++= hits
          working       ++= hits
          remaining       = misses
        }
      }

      assignIdsToNodes(roots)
    }
  }


  /**
    * Version of the adjacency assigner that works for paired UMIs stored as a single tag of
    * the form A-B where reads with A-B and B-A are related but not identical.
    *
    * @param maxMismatches the maximum number of mismatches between UMIs
    */
  class PairedUmiAssigner(maxMismatches: Int) extends AdjacencyUmiAssigner(maxMismatches) {
    /** String that is prefixed onto the UMI from the read with that maps to a lower coordinate in the genome.. */
    private[umi] val lowerReadUmiPrefix: String = ("a" * (maxMismatches+1)) + ":"

    /** String that is prefixed onto the UMI from the read with that maps to a higher coordinate in the genome.. */
    private[umi] val higherReadUmiPrefix: String = ("b" * (maxMismatches+1)) + ":"

    /** Takes a UMI of the form "A-B" and returns "B-A". */
    def reverse(umi: Umi): Umi = umi.indexOf('-') match {
      case -1 => throw new IllegalStateException(s"UMI $umi is not a paired UMI.")
      case i  =>
        val first  = umi.substring(0, i)
        val second = umi.substring(i+1, umi.length)
        second + '-' + first
    }

    /** Turns each UMI into the lexically earlier of A-B or B-A and then counts them. */
    override protected def count(umis: Seq[Umi]): Iterator[(Umi, Long)] = {
      val counter = new SimpleCounter[Umi]()
      umis.foreach { umi =>
        val reversed = reverse(umi)
        val lower = if (umi.compare(reversed) < 0) umi else reversed
        counter.count(lower)
      }

      counter.iterator
    }

    /** Returns whether or not a paired-UMI matches within mismatch tolerance. */
    override protected def matches(lhs: Umi, rhs: Umi): Boolean = super.matches(lhs, rhs) || super.matches(reverse(lhs), rhs)

    /** Takes in a map of UMI to "sentinel" UMI, and outputs a map of UMI -> Molecule ID. */
    override protected def assignIdsToNodes(roots: Seq[Node]): Map[Umi, MoleculeId] = {
      val mappings = mutable.Buffer[(Umi,MoleculeId)]()
      roots.foreach(root => {
        val id = nextId
        val ab = id + GroupReadsByUmi.TopStrandDuplex
        val ba = id + GroupReadsByUmi.BottomStrandDuplex

        mappings.append((root.umi, ab))
        mappings.append((reverse(root.umi), ba))
        root.descendants.foreach { child =>
          val childUmi    = child.umi
          val childUmiRev = reverse(child.umi)

          // If the root UMI and child UMI are more similar then presumably they originate
          // from the same pairing of UMIs, otherwise if the root UMI is more similar to the
          // reversed child UMI, they are paired with each other's inverse
          if (countMismatches(root.umi, childUmi) < countMismatches(root.umi, childUmiRev)) {
            mappings.append((childUmi, ab))
            mappings.append((childUmiRev, ba))
          }
          else {
            mappings.append((childUmi, ba))
            mappings.append((childUmiRev, ab))
          }
        }
      })

      mappings.toMap
    }

    /** Since we generate mappings for A-B and B-A whether we've seen both or not, we override
      * the main assignment method to filter out the ones that shouldn't be in the Map.
      */
    override def assign(rawUmis: Seq[Umi]): Map[Umi, MoleculeId] = {
      super.assign(rawUmis).filter { case (umi, id) => rawUmis.contains(umi) }
    }
  }
}

/**
  * Metrics produced by `GroupReadsByUmi` to describe the distribution of tag family sizes
  * observed during grouping.
  *
  * @param family_size The family size, or number of templates/read-pairs belonging to the family.
  * @param count The number of families (or source molecules) observed with `family_size` observations.
  * @param fraction The fraction of all families of all sizes that have this specific `family_size`.
  * @param fraction_gt_or_eq_family_size The fraction of all families that have `>= family_size`.
  */
case class TagFamilySizeMetric(family_size: Int,
                               count: Count,
                               var fraction: Proportion = 0,
                               var fraction_gt_or_eq_family_size: Proportion = 0) extends Metric

/** The strategies implemented by [[GroupReadsByUmi]] to identify reads from the same source molecule.*/
sealed trait Strategy extends EnumEntry {
  def newStrategy(edits: Int): UmiAssigner
}
object Strategy extends FgBioEnum[Strategy] {
  def values: IndexedSeq[Strategy] = findValues
  /** Strategy to only reads with identical UMI sequences are grouped together. */
  case object Identity extends Strategy {
    def newStrategy(edits: Int = 0): UmiAssigner = {
      require(edits == 0, "Edits should be zero when using the identity UMI assigner.")
      new IdentityUmiAssigner
    }
  }
  /** Strategy to cluster reads into groups based on mismatches between reads in clusters. */
  case object Edit extends Strategy { def newStrategy(edits: Int): UmiAssigner = new SimpleErrorUmiAssigner(edits) }
  /** Strategy based on the directed adjacency method described in [umi_tools](http://dx.doi.org/10.1101/051755)
    * that allows for errors between UMIs but only when there is a count gradient.
    */
  case object Adjacency extends Strategy { def newStrategy(edits: Int): UmiAssigner = new AdjacencyUmiAssigner(edits) }
  /** Strategy similar to the [[Adjacency]] strategy similar to adjacency but for methods that produce template with a
    * pair of UMIs such that a read with A-B is related to but not identical to a read with B-A.
    */
  case object Paired extends Strategy { def newStrategy(edits: Int): UmiAssigner = new PairedUmiAssigner(edits)}
}

@clp(group=ClpGroups.Umi, description =
  """
    |Groups reads together that appear to have come from the same original molecule. Reads
    |are grouped by template, and then templates are sorted by the 5' mapping positions of
    |the reads from the template, used from earliest mapping position to latest. Reads that
    |have the same end positions are then sub-grouped by UMI sequence.
    |
    |Accepts reads in any order (including `unsorted`) and outputs reads sorted by:
    |
    |   1. The lower genome coordinate of the two outer ends of the templates
    |   2. The sequencing library
    |   3. The assigned UMI tag
    |   4. Read Name
    |
    |Reads are aggressively filtered out so that only high quality read pairs with both ends
    |mapped are taken forward.  (Note: the `MQ` tag is required on reads with mapped mates).
    |This is done with the expectation that the next step is building consensus reads, where
    |it is undesirable to either:
    |
    |   1. Assign reads together that are really from different source molecules
    |   2. Build two groups from reads that are really from the same molecule
    |
    |Errors in mapping reads could lead to both and therefore are minimized.
    |
    |Grouping of UMIs is performed by one of three strategies:
    |
    |1. **identity**:  only reads with identical UMI sequences are grouped together. This strategy
    |                  may be useful for evaluating data, but should generally be avoided as it will
    |                  generate multiple UMI groups per original molecule in the presence of errors.
    |2. **edit**:      reads are clustered into groups such that each read within a group has at least
    |                  one other read in the group with <= edits differences and there are inter-group
    |                  pairings with <= edits differences. Effective when there are small numbers of
    |                  reads per UMI, but breaks down at very high coverage of UMIs.
    |3. **adjacency**: a version of the directed adjacency method described in [umi_tools](http://dx.doi.org/10.1101/051755)
    |                  that allows for errors between UMIs but only when there is a count gradient.
    |4. **paired**:    similar to adjacency but for methods that produce template with a pair of UMIs
    |                  such that a read with A-B is related to but not identical to a read with B-A.
    |                  Expects the pair of UMIs to be stored in a single tag, separated by a hyphen
    |                  (e.g. `ACGT-CCGG`).  The molecular IDs produced have more structure than for single
    |                  UMI strategies, and are of the form `{base}/{AB|BA}`. E.g. two UMI pairs would be
    |                  mapped as follows AAAA-GGGG -> 1/AB, GGGG-AAAA -> 1/BA.
    |
    |`edit`, `adjacency` and `paired` make use of the `--edits` parameter to control the matching of
    |non-identical UMIs.
    |
    |By default, all UMIs must be the same length. If `--min-umi-length=len` is specified then reads that have a UMI
    |shorter than `len` will be discarded, and when comparing UMIs of different lengths, the first len bases will be
    |compared, where `len` is the length of the shortest UMI. The UMI length is the number of [ACGT] bases in the UMI
    |(i.e. does not count dashes and other non-ACGT characters). This option is not implemented for reads with UMI pairs
    |(i.e. using the paired assigner).
  """
)
class GroupReadsByUmi
( @arg(flag='i', doc="The input BAM file.")              val input: PathToBam  = Io.StdIn,
  @arg(flag='o', doc="The output BAM file.")             val output: PathToBam = Io.StdOut,
  @arg(flag='f', doc="Optional output of tag family size counts.") val familySizeHistogram: Option[FilePath] = None,
  @arg(flag='t', doc="The tag containing the raw UMI.")  val rawTag: String    = "RX",
  @arg(flag='T', doc="The output tag for UMI grouping.") val assignTag: String = "MI",
  @arg(flag='m', doc="Minimum mapping quality.")         val minMapQ: Int      = 30,
  @arg(flag='n', doc="Include non-PF reads.")            val includeNonPfReads: Boolean = false,
  @arg(flag='s', doc="The UMI assignment strategy.") val strategy: Strategy,
  @arg(flag='e', doc="The allowable number of edits between UMIs.") val edits: Int = 1,
  @arg(flag='l', doc= """The minimum UMI length. If not specified then all UMIs must have the same length,
                       |otherwise discard reads with UMIs shorter than this length and allow for differing UMI lengths.
                       |""")
  val minUmiLength: Option[Int] = None
)extends FgBioTool with LazyLogging {
  import GroupReadsByUmi._

  require(this.minUmiLength.forall(_ => this.strategy != Strategy.Paired), "Paired strategy cannot be used with --min-umi-length")

  private val assigner = strategy.newStrategy(this.edits)

  type ReadPair = (SamRecord, SamRecord)

  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)
    this.familySizeHistogram.foreach(f => Io.assertCanWriteFile(f))

    val in = SamSource(input)
    val header = in.header
    val sorter = Bams.sorter(SamOrder.TemplateCoordinate, header)
    val sortProgress = ProgressLogger(logger, verb="Sorted")

    // A handful of counters for tracking reads
    var (filteredNonPf, filteredPoorAlignment, filteredNsInUmi, filterUmisTooShort, kept) = (0L, 0L, 0L, 0L, 0L)

    // Filter and sort the input BAM file
    logger.info("Filtering and sorting input.")
    in.iterator
      .filter(r => !r.secondary && !r.supplementary)
      .filter(r => (includeNonPfReads || r.pf)            || { filteredNonPf += 1; false })
      .filter(r => (r.paired && r.mapped && r.mateMapped) || { filteredPoorAlignment += 1; false })
      .filter(r => (r.refIndex == r.mateRefIndex)         || { filteredPoorAlignment += 1; false })
      .filter(r => (r.mapq >= this.minMapQ && r.get[Int](SAMTag.MQ.name()).exists(_ >= this.minMapQ)) || { filteredPoorAlignment += 1; false })
      .filter(r => !r.get[String](rawTag).exists(_.contains('N')) || { filteredNsInUmi += 1; false })
      .filter { r =>
        this.minUmiLength.forall { l =>
          r.get[String](this.rawTag).forall { umi =>
            l <= umi.toUpperCase.count(c => SequenceUtil.isUpperACGTN(c.toByte))
          }
        } || { filterUmisTooShort += 1; false}
      }
      .foreach(r => { sorter += r; kept += 1; sortProgress.record(r) })

    logger.info(f"Accepted $kept%,d reads for grouping.")
    if (filteredNonPf > 0) logger.info(f"Filtered out $filteredNonPf%,d non-PF reads.")
    logger.info(f"Filtered out $filteredPoorAlignment%,d reads that were not part of a high confidence FR mapped read pair.")
    logger.info(f"Filtered out $filteredNsInUmi%,d reads that contained one or more Ns in their UMIs.")
    this.minUmiLength.foreach { _ => logger.info(f"Filtered out $filterUmisTooShort%,d reads that contained UMIs that were too short.") }

    // Output the reads in the new ordering
    logger.info("Assigning reads to UMIs and outputting.")
    val outHeader = header.clone()
    SamOrder.TemplateCoordinate.applyTo(outHeader)
    val out = SamWriter(output, outHeader)

    val iterator = sorter.iterator.grouped(2).map { case Seq(x,y) => if (x.firstOfPair) (x, y) else (y, x) }.buffered // consume in pairs
    val tagFamilySizeCounter = new NumericCounter[Int]()

    while (iterator.hasNext) {
      // Take the next set of pairs by position and assign UMIs
      val pairs = takeNextGroup(iterator)
      pairs.foreach { case(r1, r2) =>
        assert(r1.name == r2.name && r1.firstOfPair && r2.secondOfPair, s"Reads out of order @ ${r1.id} + ${r2.id}")
      }
      assignUmiGroups(pairs)

      // Then output the records in the right order (assigned tag, read name, r1, r2)
      val pairsByAssignedTag = pairs.groupBy { case (r1,r2) => r1[String](this.assignTag) }

      pairsByAssignedTag.keys.toSeq.sortBy(id => (id.length, id)).foreach(tag => {
        pairsByAssignedTag(tag).sortBy(pair => pair._1.name).flatMap(pair => Seq(pair._1, pair._2)).foreach(rec => {
          out += rec
        })
      })

      // Count up the family sizes
      pairsByAssignedTag.values.foreach(ps => tagFamilySizeCounter.count(ps.size))
    }

    out.close()

    // Write out the family size histogram
    this.familySizeHistogram match {
      case None    => Unit
      case Some(p) =>
        val ms = tagFamilySizeCounter.toSeq.map { case (size, count) => TagFamilySizeMetric(size, count)}
        val total = ms.map(_.count.toDouble).sum
        ms.foreach(m => m.fraction = m.count / total)
        ms.tails.foreach { tail => tail.headOption.foreach(m => m.fraction_gt_or_eq_family_size = tail.map(_.fraction).sum) }
        Metric.write(p, ms)
    }
  }

  /** Consumes the next group of pairs with all matching end positions and returns them. */
  def takeNextGroup(iterator: BufferedIterator[ReadPair]) : Seq[ReadPair] = {
    val first = iterator.next()
    val firstEnds = ReadInfo(first._1)
    val buffer = ListBuffer[(SamRecord, SamRecord)]()
    while (iterator.hasNext && firstEnds == ReadInfo(iterator.head._1)) buffer += iterator.next()
    first :: buffer.toList
  }

  /**
    * Takes in pairs of reads and writes an attribute (specified by `assignTag`) denoting
    * sub-grouping into UMI groups by original molecule.
    */
  def assignUmiGroups(pairs: Seq[ReadPair]): Unit = {
    val umis    = truncateUmis(pairs.map { case (r1, r2) => umiForRead(r1, r2) })
    val rawToId = this.assigner.assign(umis)

    pairs.iterator.zip(umis.iterator).foreach { case ((r1, r2), umi) =>
      val id  = rawToId(umi)
      r1(this.assignTag) = id
      r2(this.assignTag) = id
    }
  }

  /** When a minimum UMI length is specified, truncates all the UMIs to the length of the shortest UMI.  For the paired
    * assigner, truncates the first UMI and second UMI sepeartely.*/
  private def truncateUmis(umis: Seq[Umi]): Seq[Umi] = this.minUmiLength match {
    case None => umis
    case Some(length) =>
      this.assigner match {
        case _: PairedUmiAssigner =>
          throw new IllegalStateException("Cannot used the paired umi assigner when min-umi-length is defined.")
        case _ =>
          val minLength = umis.map(_.length).min
          require(length <= minLength, s"Bug: UMI found that had shorter length than expected ($minLength < $length)")
          umis.map(_.substring(0, minLength))
      }
  }

  /**
    * Retrieves the UMI to use for a read pair.  In the case of single-umi strategies this is just
    * the upper-case UMI sequence.
    *
    * For Paired the UMI is extended to encode which "side" of the template the read came from - the earlier
    * or later read on the genome.  This is necessary to ensure that, when the two paired UMIs are the same
    * or highly similar, that the A vs. B groups are constructed correctly.
    */
  private def umiForRead(r1: SamRecord, r2: SamRecord): Umi = r1.get[String](this.rawTag) match {
    case None | Some("") => fail(s"Record '$r1' was missing the raw UMI tag '${this.rawTag}'")
    case Some(chs)       =>
      val umi = chs.toUpperCase

      this.assigner match {
        case paired: PairedUmiAssigner =>
          require(r1.refIndex == r2.refIndex, s"Mates on different references not supported: ${r1.name}")
          val pos1 = if (r1.positiveStrand) r1.unclippedStart else r1.unclippedEnd
          val pos2 = if (r2.positiveStrand) r2.unclippedStart else r2.unclippedEnd
          val r1Lower = pos1 < pos2 || pos1 == pos2 && r1.positiveStrand
          val umis = umi.split('-')
          require(umis.length == 2, s"Paired strategy used but umi did not contain 2 segments: $umi")

          if (r1Lower) paired.lowerReadUmiPrefix  + ":" + umis(0) + "-" + paired.higherReadUmiPrefix + ":" + umis(1)
          else         paired.higherReadUmiPrefix + ":" + umis(0) + "-" + paired.lowerReadUmiPrefix  + ":" + umis(1)
        case _ =>
          umi
      }
  }
}
