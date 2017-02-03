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

import java.nio.file.Paths
import java.util.concurrent.atomic.AtomicLong

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.Sequences.countMismatches
import com.fulcrumgenomics.util._
import dagr.commons.util.LazyLogging
import dagr.sopt.cmdline.ValidationException
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.util.SortingCollection

import scala.collection.mutable
import scala.collection.mutable.ListBuffer


object GroupReadsByUmi {
  /** A type representing an actual UMI that is a string of DNA sequence. */
  type Umi = String
  /** A type to represent a unique ID assigned to a molecule. */
  type MoleculeId = String

  private val ReadInfoTempAttributeName = "__GRBU_ReadInfo"

  /** A case class to represent all the information we need to order reads for duplicate marking / grouping. */
  case class ReadInfo(refIndex: Int, start1: Int, start2: Int, strand1: Boolean, strand2: Boolean, library: String) extends Ordered[ReadInfo] {
    override def compare(that: ReadInfo): Int = ReadInfo.Comparisons.map(f => f(this, that)).find(_ != 0).getOrElse(0)
  }

  object ReadInfo {
    /* The set of comparisons by which we order ReadInfo objects. */
    private val Comparisons: Seq[(ReadInfo,ReadInfo) => Int] = Seq(
      (lhs, rhs) => lhs.refIndex - rhs.refIndex,
      (lhs, rhs) => lhs.start1   - rhs.start1,
      (lhs, rhs) => lhs.start2   - rhs.start2,
      (lhs, rhs) => lhs.strand1.compare(rhs.strand1),
      (lhs, rhs) => lhs.strand2.compare(rhs.strand2),
      (lhs, rhs) => lhs.library.compare(rhs.library)
    )

    /** Looks in all the places the library name can be hiding. Returns the library name
      * if one is found, otherwise returns "unknown".
      */
    private def library(rec: SAMRecord): String = {
      val rg = rec.getReadGroup
      if (rg != null && rg.getLibrary != null) rg.getLibrary
      else "unknown"
    }

    /** Creates/retrieves a ReadEnds object from a SAMRecord and stores it in a temporary attribute for later user. */
    def apply(rec: SAMRecord) : ReadInfo = {
      val tmp = rec.getTransientAttribute(GroupReadsByUmi.ReadInfoTempAttributeName)
      if (tmp != null) {
        tmp.asInstanceOf[ReadInfo]
      }
      else {
        if (rec.getReferenceIndex != rec.getMateReferenceIndex) throw new IllegalArgumentException("Mate on different chrom.")
        val chrom   = rec.getReferenceIndex.toInt
        val recNeg  = rec.getReadNegativeStrandFlag
        val recPos  = if (recNeg) rec.getUnclippedEnd else rec.getUnclippedStart
        val mateNeg = rec.getMateNegativeStrandFlag
        val matePos = if (mateNeg) SAMUtils.getMateUnclippedEnd(rec) else SAMUtils.getMateUnclippedStart(rec)
        val lib = library(rec)

        val result = if (recPos < matePos || (recPos == matePos && !recNeg)) {
          new ReadInfo(chrom, recPos, matePos, recNeg, mateNeg, lib)
        }
        else {
          new ReadInfo(chrom, matePos, recPos, mateNeg, recNeg, lib)
        }

        rec.setTransientAttribute(GroupReadsByUmi.ReadInfoTempAttributeName, result)
        result
      }
    }
  }

    /** Comparator to sort the reads by the earlier of the unclipped 5' position of the
    * first or second read.  Groups reads in a convenient way for duplicate-marking and
    * UMI assignment, and ensures that both ends of a read come together.
    *
    * Private because it has some serious restrictions!  Only allows:
    *   - Paired reads only
    *   - Mapped reads with mapped mates only
    *   - No secondary or supplementary reads
    *   - Read and mate _must_ be mapped to the same chromosome
    * If any of these conditions are violated it will go badly!
    */
  private[umi] class EarlierReadComparator extends SAMRecordComparator with Ordering[SAMRecord] {
    override def fileOrderCompare(lhs: SAMRecord, rhs: SAMRecord): Int = compare(lhs, rhs)

    /** Compares two reads for sort order. */
    override def compare(lhs: SAMRecord, rhs: SAMRecord): Int = {
      // Do some asserting!
      assertValidRead(lhs)
      assertValidRead(rhs)
      val l = ReadInfo(lhs)
      val r = ReadInfo(rhs)
      val result = l.compare(r)

      if (result == 0) lhs.getReadName.compareTo(rhs.getReadName)
      else result
    }

    /** Asserts that we didn't get reads we are unable to sort. */
    final def assertValidRead(rec: SAMRecord): Unit = {
      assert(rec.getReadPairedFlag,    "Unpaired read: " + rec.getReadName)
      assert(!rec.getReadUnmappedFlag, "Unmapped read: " + rec.getReadName)
      assert(!rec.getMateUnmappedFlag, "Read w/unmapped mate: " + rec.getReadName)
      assert(rec.getReferenceIndex == rec.getMateReferenceIndex, "Read w/mate on different chr: " + rec.getReadName)
      assert(SAMUtils.getMateCigarString(rec) != null, "Read w/o Mate Cigar tag: " + rec.getReadName)
      assert(rec.getAttribute("MQ") != null, "Read w/o Mate MQ tag: " + rec.getReadName)
      assert(!rec.isSecondaryOrSupplementary, "Secondary or supplementary read: " + rec.getReadName)
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
    }

    /** Returns the count of each raw UMI that was observed. */
    protected def count(umis: Seq[Umi]): Iterator[(Umi,Long)] = SimpleCounter(umis).iterator

    /** Returns the number of differences between a pair of UMIs. */
    protected def differences(lhs: Umi, rhs: Umi): Int = Sequences.countMismatches(lhs, rhs)

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
      // Make a list of counts of all UMIs in order from most to least abundant
      val nodes = count(rawUmis).map{case(umi,count) => new Node(umi, count)}.toBuffer.sortBy((n:Node) => -n.count)
      val roots = mutable.Buffer[Node]()

      /** Function that takes a root or starting node and finds all children in the "others"
        * that are within a fixed edit mismatch distance (specified by this.maxMismatches) and
        * where the child has a count that is more plausibly explained by being an error from
        * the root than by being a different UMI. The relationship used in UMI tools is:
        *   count(root) >= 2 * count(child) - 1
        *
        * This allows a root with a single observation to have a child with a single observation,
        * but as counts of the root go up they must be >= approximately twice the child count.
        */
      def addChildren(root: Node, others: mutable.Buffer[Node]) : Unit = {
        val xs = others.filter(other => root.count >= 2 * other.count - 1 && differences(root.umi, other.umi) <= maxMismatches )
        root.children ++= xs
        others --= xs
        root.children.foreach(r => addChildren(r, others))
      }

      // Now build one or more graphs starting with the most abundant remaining umi
      while (nodes.nonEmpty) {
        val root = nodes.remove(0)
        addChildren(root, nodes)
        roots += root
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
    /** Takes a UMI of the form "A-B" and returns "B-A". */
    def reverse(umi: Umi): Umi = umi.indexOf('-') match {
      case -1 => throw new IllegalStateException(s"UMI ${umi} is not a paired UMI.")
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

    /** Returns the differences between a pair of UMIs. */
    override protected def differences(lhs: Umi, rhs: Umi): Int = Math.min(countMismatches(lhs, rhs), countMismatches(reverse(lhs), rhs))

    /** Takes in a map of UMI to "sentinel" UMI, and outputs a map of UMI -> Molecule ID. */
    override protected def assignIdsToNodes(roots: Seq[Node]): Map[Umi, MoleculeId] = {
      val mappings = mutable.Buffer[(Umi,MoleculeId)]()
      roots.foreach(root => {
        val id = nextId
        val ab = id + "/A"
        val ba = id + "/B"

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

/** Trivial metrics class for representing a tag family size histogram entry. */
case class TagFamilySizeMetric(family_size: Int,
                               count: Long,
                               var fraction: Double = 0,
                               var fraction_gt_or_eq_family_size: Double = 0) extends Metric

@clp(group=ClpGroups.Umi, description =
  """
    |Groups reads together that appear to have come from the same original molecule. Reads
    |are grouped by template, and then templates are sorted by the 5' mapping positions of
    |the reads from the template, used from earliest mapping position to latest. Reads that
    |have the same end positions are then sub-grouped by UMI sequence.
    |
    |Accepts reads in any order (including unsorted) and outputs reads sorted by:
    |   1. The lower genome coordinate of the two outer ends of the templates
    |   2. The sequencing library
    |   3. The assigned UMI tag
    |   4. Read Name
    |
    |Reads are aggressively filtered out so that only high quality read pairs with both ends
    |mapped are taken forward.  (Note: the MQ tag is required on reads with mapped mates).
    |This is done with the expectation that the next step is building consensus reads, where
    |it is undesirable to either:
    |   1. Assign reads together that are really from different source molecules
    |   2. Build two groups from reads that are really from the same molecule
    |Errors in mapping reads could lead to both and therefore are minimized.
    |
    |Grouping of UMIs is performed by one of three strategies:
    |   1. identity:  only reads with identical UMI sequences are grouped together. This strategy
    |                 may be useful for evaluating data, but should generally be avoided as it will
    |                 generate multiple UMI groups per original molecule in the presence of errors.
    |   2. edit:      reads are clustered into groups such that each read within a group has at least
    |                 one other read in the group with <= edits differences and there are inter-group
    |                 pairings with <= edits differences. Effective when there are small numbers of
    |                 reads per UMI, but breaks down at very high coverage of UMIs.
    |   3. adjacency: a version of the directed adjacency method described in umi_tools
    |                 (ref: http://dx.doi.org/10.1101/051755).
    |   4. paired:    similar to adjacency but for methods that produce template with a pair of UMIs
    |                 such that a read with A-B is related to but not identical to a read with B-A.
    |                 Expects the pair of UMIs to be stored in a single tag, separated by a hyphen
    |                 (e.g. ACGT-CCGG).  The molecular IDs produced have more structure than for single
    |                 UMI strategies, and are of the form {base}/{AB|BA}.  E.g. two UMI pairs would be
    |                 mapped as follows AAAA-GGGG -> 1/AB, GGGG-AAAA -> 1/BA.
    |
    |'edit', 'adjacency' and 'paired' make use of the 'edits' parameter to control the matching of
    |non-identical UMIs.
  """
)
class GroupReadsByUmi
( @arg(flag="i", doc="The input BAM file.")              val input: PathToBam  = Io.StdIn,
  @arg(flag="o", doc="The output BAM file.")             val output: PathToBam = Io.StdOut,
  @arg(flag="f", doc="Optional output of tag family size counts.") val familySizeHistogram: Option[FilePath] = None,
  @arg(flag="t", doc="The tag containing the raw UMI.")  val rawTag: String    = "RX",
  @arg(flag="T", doc="The output tag for UMI grouping.") val assignTag: String = "MI",
  @arg(flag="m", doc="Minimum mapping quality.")         val minMapQ: Int      = 30,
  @arg(flag="n", doc="Include non-PF reads.")            val includeNonPfReads: Boolean = false,
  @arg(flag="s", doc="The UMI assignment strategy; one of 'identity', 'edit', 'adjacency' or 'paired'.") val strategy: String,
  @arg(flag="e", doc="The allowable number of edits between UMIs.") val edits: Int = 1,
  @arg(          doc="Temporary directory for sorting.") val tmpDir: DirPath = Paths.get(System.getProperty("java.io.tmpdir"))
)extends FgBioTool with LazyLogging {
  import GroupReadsByUmi._

  private val assigner = strategy match {
    case "identity"  => new IdentityUmiAssigner
    case "edit"      => new SimpleErrorUmiAssigner(edits)
    case "adjacency" => new AdjacencyUmiAssigner(edits)
    case "paired"    => new PairedUmiAssigner(edits)
    case other       => throw new ValidationException(s"Unknown strategy: $other")
  }

  type ReadPair = (SAMRecord, SAMRecord)

  override def execute(): Unit = {
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)
    this.familySizeHistogram.foreach(f => Io.assertCanWriteFile(f))

    val in = SamReaderFactory.make().open(input.toFile)
    val header = in.getFileHeader
    val sorter = SortingCollection.newInstance(classOf[SAMRecord], new BAMRecordCodec(header), new EarlierReadComparator, 100000, tmpDir.toFile)
    val sortProgress = new ProgressLogger(logger, verb="Sorted")

    // Filter and sort the input BAM file
    logger.info("Filtering and sorting input.")
    in.iterator
      .filter(r => includeNonPfReads || !r.getReadFailsVendorQualityCheckFlag)
      .filter(r => !r.isSecondaryOrSupplementary)
      .filter(r => r.getReadPairedFlag)
      .filter(r => !r.getReadUnmappedFlag && !r.getMateUnmappedFlag)
      .filter(r => r.getReferenceIndex == r.getMateReferenceIndex)
      .filter(r => r.getMappingQuality >= this.minMapQ && Option(r.getIntegerAttribute(SAMTag.MQ.name())).exists(_ >= this.minMapQ))
      .foreach(r => { sorter.add(r); sortProgress.record(r) })

    // Output the reads in the new ordering
    logger.info("Assigning reads to UMIs and outputting.")
    val outHeader = header.clone()
    outHeader.setSortOrder(SortOrder.unsorted)
    outHeader.setGroupOrder(GroupOrder.query)
    val out = new SAMFileWriterFactory().makeWriter(outHeader, true, output.toFile, null)
    val outProgress = new ProgressLogger(logger, verb="Output")

    val iterator = sorter.iterator().grouped(2).map(i => (i(0), i(1))).buffered // consume in pairs
    val tagFamilySizeCounter = new NumericCounter[Int]()

    while (iterator.hasNext) {
      // Take the next set of pairs by position and assign UMIs
      val pairs = takeNextGroup(iterator)
      pairs.foreach { case(r1, r2) => assert(r1.getReadName == r2.getReadName, "Reads out of order @ " + r1.getReadName) }
      assignUmiGroups(pairs)

      // Then output the records in the right order (assigned tag, read name, r1, r2)
      val pairsByAssignedTag = pairs.groupBy(pair => pair._1.getStringAttribute(this.assignTag))

      pairsByAssignedTag.keys.toSeq.sorted.foreach(tag => {
        pairsByAssignedTag(tag).sortBy(pair => pair._1.getReadName).flatMap(pair => Seq(pair._1, pair._2)).foreach(rec => {
          out.addAlignment(rec)
          outProgress.record(rec)
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
        Metric.write(ms, p)
    }
  }

  /** Consumes the next group of pairs with all matching end positions and returns them. */
  def takeNextGroup(iterator: BufferedIterator[ReadPair]) : Seq[ReadPair] = {
    val first = iterator.next()
    val firstEnds = ReadInfo(first._1)
    val buffer = ListBuffer[(SAMRecord, SAMRecord)]()
    while (iterator.hasNext && firstEnds == ReadInfo(iterator.head._1)) buffer += iterator.next()
    first :: buffer.toList
  }

  /**
    * Takes in pairs of reads and writes an attribute (specified by `assignTag`) denoting
    * sub-grouping into UMI groups by original molecule.
    */
  def assignUmiGroups(pairs: Seq[ReadPair]): Unit = {
    val rawUmis = pairs.map(_._1).map(getRawUmi)
    val rawToId = this.assigner.assign(rawUmis)

    pairs.foreach(pair => Seq(pair._1, pair._2).foreach(rec => {
      val raw = getRawUmi(rec)
      val id  = rawToId(raw)
      rec.setAttribute(this.assignTag, id)
    }))
  }

  private def getRawUmi(rec: SAMRecord): String = rec.getStringAttribute(this.rawTag) match {
    case null | "" => fail(s"Record '$rec' was missing the raw UMI tag '${this.rawTag}'")
    case value => value.toUpperCase
  }
}
