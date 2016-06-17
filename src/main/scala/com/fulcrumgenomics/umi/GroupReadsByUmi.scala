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

import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.{Io, ProgressLogger, Sequences, SimpleCounter}
import com.fulcrumgenomics.util.Sequences.countMismatches
import dagr.commons.CommonsDef.{DirPath, PathToBam}
import dagr.commons.util.LazyLogging
import dagr.sopt.cmdline.ValidationException
import dagr.sopt.{arg, clp}
import htsjdk.samtools.SAMFileHeader.{GroupOrder, SortOrder}
import htsjdk.samtools._
import htsjdk.samtools.util.SortingCollection

import scala.collection.mutable
import scala.collection.mutable.ListBuffer


object GroupReadsByUmi {
  private val SortPositionAttribute        = "__GRs"
  private val LibraryAttribute             = "__GRlb"
  private val InsertEndsAndStrandAttribute = "__Insert"

  /** Trait that can be implemented to provide a UMI assignment strategy. */
  private[umi] sealed trait UmiAssigner {
    def assign(rawUmis: Seq[String]) : Map[String, String]
  }

  /** Assigns UMIs based solely on sequence identity. */
  private[umi] class IdentityUmiAssigner extends UmiAssigner {
    override def assign(rawUmis: Seq[String]): Map[String, String] = rawUmis.map(umi => umi -> umi.toUpperCase).toMap
  }

  /** Assigns UMIs based solely on sequence identity. */
  private[umi] class SimpleErrorUmiAssigner(val maxMismatches: Int) extends UmiAssigner {
    override def assign(rawUmis: Seq[String]): Map[String, String] = {
      type UmiSet = mutable.SortedSet[String]
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

      umiSets.flatMap(set => set.map(umi => umi -> set.head)).toMap
    }
  }

  /**
    * Class that implements the directed adjacency graph method from umi_tools.
    * See: https://github.com/CGATOxford/UMI-tools
    */
  private[umi] class AdjacencyUmiAssigner(val maxMismatches: Int) extends UmiAssigner {
    /** Represents a node in the adjacency graph; equality is just by UMI sequence. */
    private class Node(val umi: String, val count: Long, val children: mutable.Buffer[Node] = mutable.Buffer()) {
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

    override def assign(rawUmis: Seq[String]): Map[String, String] = {
      // Make a list of counts of all UMIs in order from most to least abundant
      val nodes = SimpleCounter(rawUmis).map{case(umi,count) => new Node(umi, count)}.toBuffer.sortBy((n:Node) => -n.count)
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
        val xs = others.filter(other => countMismatches(root.umi, other.umi) <= maxMismatches && root.count >= 2 * other.count - 1)
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

      // And finally make the output mapping
      val mappings = mutable.Buffer[(String,String)]()
      roots.foreach(root => {
        mappings.append((root.umi, root.umi))
        root.descendants.foreach(child => mappings.append((child.umi, root.umi)))
      })

      mappings.toMap
    }
  }

  /** Looks in all the places the library name can be hiding. Returns the library name
    * if one is found, otherwise returns "unknown".
    */
  final def library(rec: SAMRecord): String = {
    var lib: String = rec.getTransientAttribute(LibraryAttribute).asInstanceOf[String]
    if (lib == null) {
      val rg = rec.getReadGroup
      if (rg != null && rg.getLibrary != null) lib = rg.getLibrary
      else lib = rec.getStringAttribute("LB")

      if (lib == null) lib = "unknown"
      rec.setTransientAttribute(LibraryAttribute, lib)
    }

    lib
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
  private[umi] class EarlierReadComparator extends SAMRecordComparator {
    override def fileOrderCompare(lhs: SAMRecord, rhs: SAMRecord): Int = compare(lhs, rhs)

    /** A set of functions that extract a value from a read, used in order to compare reads. */
    private val comparisons: Seq[(SAMRecord,SAMRecord) => Int] = Seq(
      (lhs, rhs) => lhs.getReferenceIndex - rhs.getReferenceIndex,
      (lhs, rhs) => sortPosition(lhs) - sortPosition(rhs),
      (lhs, rhs) => library(lhs).compareTo(library(rhs)),
      (lhs, rhs) => lhs.getReadName.compareTo(rhs.getReadName),
      (lhs, rhs) => if (lhs.getFirstOfPairFlag) -1 else 1
    )

    /** Compares two reads for sort order. */
    override def compare(lhs: SAMRecord, rhs: SAMRecord): Int = {
      // Do some asserting!
      assertValidRead(lhs)
      assertValidRead(rhs)
      comparisons.iterator.map(f => f(lhs, rhs)).find(_ != 0).getOrElse(0)
    }

    /** Asserts that we didn't get reads we are unable to sort. */
    final def assertValidRead(rec: SAMRecord) = {
      assert(rec.getReadPairedFlag,    "Unpaired read: " + rec.getReadName)
      assert(!rec.getReadUnmappedFlag, "Unmapped read: " + rec.getReadName)
      assert(!rec.getMateUnmappedFlag, "Read w/unmapped mate: " + rec.getReadName)
      assert(rec.getReferenceIndex == rec.getMateReferenceIndex, "Read w/mate on different chr: " + rec.getReadName)
      assert(SAMUtils.getMateCigarString(rec) != null, "Read w/o Mate Cigar tag: " + rec.getReadName)
      assert(rec.getAttribute("MQ") != null, "Read w/o Mate MQ tag: " + rec.getReadName)
      assert(!rec.isSecondaryOrSupplementary, "Secondary or supplementary read: " + rec.getReadName)
    }

    /** Retrieves the lower of the rec's and mate's 5' unclipped position for sorting. */
    final def sortPosition(rec: SAMRecord): Int = {
      val tmp = rec.getTransientAttribute(SortPositionAttribute)
      if (tmp != null) {
        tmp.asInstanceOf[Int]
      }
      else {
        val p1 = if (rec.getReadNegativeStrandFlag) rec.getUnclippedEnd else rec.getUnclippedStart
        val p2 = if (rec.getMateNegativeStrandFlag) SAMUtils.getMateUnclippedEnd(rec) else SAMUtils.getMateUnclippedStart(rec)
        val pick = Math.min(p1, p2)
        rec.setTransientAttribute(SortPositionAttribute, pick)
        pick
      }
    }
  }
}

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
    |mapped are taken forward.  This is done with the expectation that the next step is building
    |consensus reads, where it is undesirable to either:
    |   1. Assign reads together that are really from different source molecules
    |   2. Build two groups from reads that are really from the same molecule
    |Errors in mapping reads could lead to both and therefore are minimized.
    |
    |Grouping of UMIs is performed by one of three strategies:
    |   1. identity:  only reads with identical UMI sequences are grouped together. This strategy
    |                 may be usful for evaluating data, but should generally be avoided as it will
    |                 generate multiple UMI groups per original molecule in the presence of errors.
    |   2. edit:      reads are clustered into groups such that each read within a group has at least
    |                 one other read in the group with <= edits differences and there are inter-group
    |                 pairings with <= edits differences. Effective when there are small numbers of
    |                 reads per UMI, but breaks down at very high coverage of UMIs.
    |   3. adjacency: a version of the directed adjacency method described in umi_tools
    |                 (ref: http://dx.doi.org/10.1101/051755).
    |
    |Both 'edit' and 'adjacency' make use of the 'edits' parameter to control the matching of
    |non-identical UMIs.
  """
)
class GroupReadsByUmi
( @arg(flag="i", doc="The input BAM file.")              input: PathToBam  = Io.StdIn,
  @arg(flag="o", doc="The output BAM file.")             output: PathToBam = Io.StdOut,
  @arg(flag="t", doc="The tag containing the raw UMI.")  rawTag: String    = "RX",
  @arg(flag="T", doc="The output tag for UMI grouping.") assignTag: String = "MI",
  @arg(flag="m", doc="Minimum mapping quality.")         minMapQ: Int      = 30,
  @arg(flag="s", doc="The UMI assignment strategy; one of 'identity', 'edit', 'adjacency'.") strategy: String,
  @arg(flag="e", doc="The allowable number of edits between UMIs.") edits: Int = 1,
  @arg(          doc="Temporary directory for sorting.") tmpDir: DirPath = Paths.get(System.getProperty("java.io.tmpdir"))
)extends FgBioTool with LazyLogging {
  import GroupReadsByUmi._

  private val assigner = strategy match {
    case "identity"  => new IdentityUmiAssigner
    case "edit"      => new SimpleErrorUmiAssigner(edits)
    case "adjacency" => new AdjacencyUmiAssigner(edits)
    case other       => throw new ValidationException(s"Unknown strategy: $other")
  }

  private val counter = new AtomicLong()
  type ReadPair = (SAMRecord, SAMRecord)

  override def execute(): Unit = {
    import scala.collection.JavaConversions.asScalaIterator
    Io.assertReadable(input)
    Io.assertCanWriteFile(output)

    val in = SamReaderFactory.make().open(input.toFile)
    val header = in.getFileHeader
    val sorter = SortingCollection.newInstance(classOf[SAMRecord], new BAMRecordCodec(header), new EarlierReadComparator, 100000, tmpDir.toFile)
    val sortProgress = new ProgressLogger(logger, verb="Sorted")

    // Filter and sort the input BAM file
    logger.info("Filtering and sorting input.")
    in.iterator
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
    }

    out.close()
  }

  /** Consumes the next group of pairs with all matching end positions and returns them. */
  def takeNextGroup(iterator: BufferedIterator[ReadPair]) : Seq[ReadPair] = {
    val first = iterator.next()
    val firstEnds = extractEnds(first._1)
    val buffer = ListBuffer[(SAMRecord, SAMRecord)]()
    while (iterator.hasNext &&
           library(first._1) == library(iterator.head._1) &&
           firstEnds == extractEnds(iterator.head._1)) buffer += iterator.next()
    first :: buffer.toList
  }

  /**
    * Takes in pairs of reads and writes an attribute (specified by `assignTag`) denoting
    * sub-grouping into UMI groups by original molecule.
    */
  def assignUmiGroups(pairs: Seq[ReadPair]): Unit = {
    val rawUmis       = pairs.map(_._1).map(_.getStringAttribute(this.rawTag).toUpperCase)
    val rawToAssigned = this.assigner.assign(rawUmis)
    val assignedToId  = rawToAssigned.values.toSet[String].map(_ -> counter.getAndIncrement().toString).toMap

    pairs.foreach(pair => Seq(pair._1, pair._2).foreach(rec => {
      val raw = rec.getStringAttribute(this.rawTag).toUpperCase
      val assigned = rawToAssigned(raw)
      val id = assignedToId(assigned)
      rec.setAttribute(this.assignTag, id)
    }))
  }

  /** Extracts a tuple of:
    *    (RefIndex, Lower 5' End, Higher 5' End, Strand of Lower, Strand of Higher)
    * that is appropriate to group reads by for duplicate marking/UMI grouping.
    */
  def extractEnds(rec: SAMRecord) : (Int,Int,Int,Boolean,Boolean) = {
    val tmp = rec.getTransientAttribute(GroupReadsByUmi.InsertEndsAndStrandAttribute)
    if (tmp != null) {
      tmp.asInstanceOf[(Int,Int,Int,Boolean,Boolean)]
    }
    else {
      val chrom   = rec.getReferenceIndex.toInt
      val recNeg  = rec.getReadNegativeStrandFlag
      val recPos  = if (recNeg) rec.getUnclippedEnd else rec.getUnclippedStart
      val mateNeg = rec.getMateNegativeStrandFlag
      val matePos = if (mateNeg) SAMUtils.getMateUnclippedEnd(rec) else SAMUtils.getMateUnclippedStart(rec)

      val result = if (recPos < matePos || (recPos == matePos && !recNeg)) {
        (chrom, recPos, matePos, recNeg, mateNeg)
      }
      else {
        (chrom, matePos, recPos, mateNeg, recNeg)
      }

      rec.setTransientAttribute(GroupReadsByUmi.InsertEndsAndStrandAttribute, result)
      result
    }
  }
}
