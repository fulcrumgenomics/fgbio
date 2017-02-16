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
 */

package com.fulcrumgenomics.util

import java.util.Random
import java.util.regex.Pattern

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.util.PickLongIndices.Index
import dagr.commons.util.LazyLogging
import dagr.sopt.{arg, clp}
import htsjdk.samtools.util.SequenceUtil

import scala.annotation.tailrec
import scala.collection.mutable
import scala.collection.mutable.ListBuffer


object PickLongIndices {
  type Index = Array[Byte]
}

/** A set of indices that maintains things as both Strings and byte[]s. */
private class IndexSet extends Iterable[Index] {
  private val arrays = ListBuffer[Array[Byte]]()
  private val set    = mutable.HashSet[String]()

  override def iterator: Iterator[Array[Byte]] = this.arrays.iterator

  /** Add an index to the set. */
  def +=(index : Array[Byte]): Unit = {
    this.arrays += index
    this.set.add(new String(index))
  }

  /** Checks to see whether an index is contained in the set. */
  def contains(index: Array[Byte]): Boolean = contains(new String(index))

  /** Checks to see whether an index is contained in the set. */
  def contains(index: String): Boolean = this.set.contains(index)
}

/** A Metric class for information about each picked index. */
case class IndexMetric(index: String,
                       source: String,
                       min_mismatches: Int,
                       indices_at_min_mismatches: Int,
                       gc: Double,
                       longest_homopolymer: Int,
                       worst_structure_seq: Option[String],
                       worst_structure_dbn: Option[String],
                       worst_structure_delta_g: Option[Double]) extends Metric {

  override protected def formatValues(value: Any): String = value match {
    case None    => ""
    case Some(x) => super.formatValues(x)
    case x       => super.formatValues(x)
  }
}

@clp(group=ClpGroups.Utilities, description =
    """
      |Picks a set of molecular indices that have at least a given number of mismatches between
      |them. Whereas PickIlluminaIndices attempts to pick a near-optimal set of indices,
      |PickLongIndices implements a significantly more efficient method based on generation of
      |random indices that can generate a large set of satisfactory indices in a small amount of
      |time and memory even for index lengths >> 10bp.
      |
      |Many options exist for controlling aspects of the indices picked, including length, edit
      |distance (mismatches only), gc range, homopolymer content, secondary structure etc.
      |
      |Secondary structure is predicted using ViennaRNA's RNAfold in DNA mode. To enable structure
      |checking both the --vienna-rna-dir and --adapters must be specified.  Adapters must be
      |strings of A, C, G, and T with a single block of Ns (e.g. ACGTNNNN or ACNNNNGT).  At runtime
      |the Ns are replaced with indices, and deltaG of the index-containing sequence is calculated.
      |
      |The number of indices requested may not be possible to produce given other constraints.
      |When this is the case the tool will output as many indices as possible, though less than
      |the requested number.  In such cases it may be useful to try different values of --attempts.
      |This parameter controls how many attempts are made to find the next valid index before
      |quitting and outputting the accumulated indices.  Higher values will yield incrementally more
      |indices but require significantly longer runtimes.
      |
      |A file of existing indices may be provided. Existing indices must be of the same length as
      |requested indices and composed of A, C, G and Ts, but are subject to no other constraints.
      |Index picking will then built a set comprised of the existing indices, and new indices which
      |satisfy all constraints.  Existing indices are included in the generated output file.
    """)
class PickLongIndices
(
  @arg(flag="l", doc="The length of each index sequence.")                             val length: Int = 8,
  @arg(flag="n", doc="The number of indices desired.")                                 val numberOfIndices: Int,
  @arg(flag="e", doc="The minimum edit distance between two indices in the set.")      val editDistance: Int = 3,
  @arg(flag="o", doc="File to write indices to.")                                      val output: FilePath,
  @arg(          doc="Allow indices that are lexical reverses of one another")         val allowReverses: Boolean = false,
  @arg(          doc="Allow indices that are reverse complements of one another")      val allowReverseComplements: Boolean = false,
  @arg(          doc="Allow indices that are palindromic (index == revcomp(index)).")  val allowPalindromes: Boolean = false,
  @arg(          doc="Reject indices with a homopolymer of greater than this length.") val maxHomopolymer: Int = 2,
  @arg(flag="g", doc="The minimum GC fraction for an index to be accepted.")           val minGc: Double = 0.2,
  @arg(flag="G", doc="The maximum GC fraction for an index to be accepted.")           val maxGc: Double = 0.8,
  @arg(          doc="File of existing index sequences to integrate, one per line.")   val existing: Option[FilePath] = None,
  @arg(flag="s", doc="Random seed value.")                                             val seed: Int = 1,
  @arg(flag="a", doc="Attempts to pick the next index before quitting.")               val attempts: Int = 1e5.toInt,
  @arg(          doc="The installation directory for ViennaRNA.")                      val viennaRnaDir: Option[DirPath] = None,
  @arg(flag="t", doc="The temperature at which to predict secondary structure.")       val temperature: Double = 25d,
  @arg(          doc="The lowest acceptable secondary structure deltaG.")              val minDeltaG: Double = -10,
  @arg(          doc="Adapter sequence(s) into which indices will be inserted.", minElements=0) val adapters: Seq[String] = Seq(),
  @arg(          doc="Any index sequence that appears in an avoid sequence or its reverse complement will be discarded.")
  val avoidSequence: Seq[String] = IlluminaAdapters.all.flatMap(_.both)
) extends FgBioTool with LazyLogging {
  // A few "constants"
  private val Bases  = Array[Byte]('A', 'C', 'G', 'T')
  private val random = new Random(seed)

  // Input checking
  existing.foreach(Io.assertReadable)
  Io.assertCanWriteFile(output)
  if (length < 1)                invalid("Length must be a positive integer.")
  if (numberOfIndices < 1)       invalid("Number of indices to pick must be a positive integer.")
  if (maxHomopolymer < 0)        invalid("Max homopolymer length must be >= 0.")
  if (minGc >= maxGc)            invalid("Min GC must be greater than max GC.")
  if (editDistance < 0)          invalid("Edit distance must be >= 0.")
  if (editDistance > length)     invalid("Edit distance must be <= length.")
  if (viennaRnaDir.isDefined && adapters.isEmpty) invalid("Specifying Vienna RNA dir has no effect if no adapters are provided.")
  if (viennaRnaDir.isEmpty && adapters.nonEmpty)  invalid("Specifying adapters has no effect if Vienna RNA dir is not provided.")

  /** Parse and check any set of existing indices. */
  private val existingIndices = existing match {
    case None    => Seq.empty[Index]
    case Some(p) => Io.readLines(p).map(_.trim.toUpperCase).filter(_.nonEmpty).map(_.getBytes).toSeq
  }

  if (existingIndices.exists(_.length != length)) invalid(s"Existing indices are not all of length ${length}")
  if (existingIndices.exists(idx => idx.exists(b => !Bases.contains(b)))) invalid("One or more existing indices contains a non-[ACGT] character.")

  /** The set of KMERs from the avoid sequences, that should not be selected. */
  private val avoidKmers = avoidSequence.flatMap { seq =>
    (0 to seq.length - length).flatMap { i =>
      val sub = seq.substring(i, i + length).toUpperCase
      Seq(sub, SequenceUtil.reverseComplement(sub))
    }
  }.toSet

  /* An Option[DnaFoldPredictor] to use to predict secondary structure. */
  private val dnaFoldPredictor = viennaRnaDir.map(dir => new DnaFoldPredictor(dir.toFile, temperature))
  private[util] val adapterRegex = Pattern.compile("^[ACGT]*N+[ACGT]*$")
  adapters.map(_.toUpperCase).filterNot(adapterRegex.matcher(_).matches()) match {
    case Seq() => Unit
    case xs    => invalid(s"Adapters must be all [ACGT] with a single contiguous block of Ns where the index goes: ${xs}")
  }

  /** The main method! */
  override def execute(): Unit = {
    val picks = pickIndices
    logger.info(s"Picked ${picks.size} indices.")
    writeOutput(picks, output)
  }

  /** Picks a set of indices. */
  private[util] def pickIndices: IndexSet = {
    val indices = new IndexSet
    existingIndices.foreach(i => indices += i)

    var ok = true
    var startTime = System.currentTimeMillis()
    while (ok && indices.size < numberOfIndices) {
      val pick = nextIndex(indices, this.length, this.editDistance, this.attempts)
      pick.foreach(p => indices += p)
      ok = pick.nonEmpty

      if (System.currentTimeMillis() - startTime > 30000) {
        logger.info(s"Picked ${indices.size} indices so far.")
        startTime = System.currentTimeMillis()
      }
    }

    indices
  }

  /**
    * Attempts to generate the next index that can be added to the set without violating any
    * of the constraints.  Will evaluate 'attempts' indices which would be acceptable independently
    * to see if they can be added.
    *
    * @return Some(Index) if an index can be found within the number of attempts, otherwise None
    */
  private def nextIndex(picks: IndexSet, length: Int, minEdits: Int, attempts: Int): Option[Index] = {
    var pick: Option[Index] = None
    var attempt = 0
    forloop (0)(_ < attempts && pick.isEmpty)(_ + 1) { attempt =>
      val index = randomAcceptableIndex(length)
      val string = new String(index)

      // structure checking (worstDeltaG) is done here instead of in randomAcceptableIndex because,
      // even though it's a per-index property, it's much slower than all the other checks, and
      // when making many attempts to find a suitable next index, the time spent structure checking
      // all the indices we're going to discard for too-few-edits anyway explodes.
      if ( !picks.contains(string) &&
           (allowReverses || !picks.contains(string.reverse)) &&
           (allowReverseComplements || !picks.contains(SequenceUtil.reverseComplement(string))) &&
           picks.forall(p => mismatches(index, p) >= minEdits) &&
           worstStructure(index).forall(_.deltaG >= minDeltaG)) {
        pick = Some(index)
      }
    }

    pick
  }

  /** Generates the next random index which is acceptable as a standalone index, prior to any secondary structure checking. */
  @tailrec private def randomAcceptableIndex(length: Int): Index = {
    val index  = randomIndex(length)
    val string = new String(index)
    val gc     = SequenceUtil.calculateGc(index)

    if (gc >= minGc && gc <= maxGc &&
        homopolymerLengthOk(index, maxHomopolymer) &&
        !avoidKmers.contains(string) &&
        (allowPalindromes || !isPalindrome(index))
    ) {
      index
    }
    else {
      randomAcceptableIndex(length)
    }
  }

  /** Returns true if a sequence is a palindrome, otherwise false. */
  private[util] def isPalindrome(seq: Index): Boolean = {
    var i = 0
    var j = seq.length-1
    var palindrome = true
    while (palindrome && i < j) {
      palindrome = SequenceUtil.basesEqual(seq(i), SequenceUtil.complement(seq(j)))
      i += 1
      j -= 1
    }
    palindrome
  }

  /** If we have structure prediction parameters, returns Some(min deltaG) else None. */
  private[util] def worstStructure(index: Index): Option[DnaFoldPrediction] = {
    val sIndex = new String(index)
    (this.dnaFoldPredictor, adapters) match {
      case (None, _)     => None
      case (_, Seq())    => None
      case (Some(f), as) => adapters.map(a => f.predict(a.replaceAll("N+", sIndex))).sortBy(_.deltaG()).headOption
    }
  }

  /** Generates a random index sequence of the desired length. */
  private def randomIndex(length: Int): Index = {
    val index = new Array[Byte](length)
    forloop (0, length) { i => index(i) = Bases(random.nextInt(4)) }
    index
  }

  /** Checks to see if the index sequence contains any homopolymers longer than longest. */
  private[util] def homopolymerLengthOk(index: Index, longest: Int): Boolean = {
    var end = index.length - longest
    var ok = true

    // Check for a homopolymer starting at each base up to longest+1 from the end of the index
    forloop (0)(ok && _ < end)(_ + 1) { i =>
      val base = index(i)
      var inHomoPolymer = true
      forloop (i+1)(inHomoPolymer && _ < i+longest+1)(_ + 1) { j =>
        if (base != index(j)) inHomoPolymer = false
      }

      ok = !inHomoPolymer
    }

    ok
  }

  /** Counts mismatches in two arrays of equal length, all upper-case non-ambiguous, bases. */
  private[util] def mismatches(lhs: Index, rhs: Index): Int = {
    var n = 0
    forloop (0, lhs.length) { i => if (lhs(i) != rhs(i)) n += 1 }
    n
  }

  /** Writes information about each index to a file. */
  private def writeOutput(indices: IndexSet, output: FilePath): Unit = {
    // Go through and compute pair-wise distance between all the indices
    case class AnnotatedIndex(index: Index, existing: Boolean, distances: NumericCounter[Int] = new NumericCounter[Int])
    val annotated = indices.map(index => AnnotatedIndex(index, existing=existingIndices.contains(index)))
      .toList.sortBy((a: AnnotatedIndex) => (!a.existing, new String(a.index)))

    annotated.tails.foreach {
      case x :: ys => ys.foreach { y =>
        val n = mismatches(x.index, y.index)
        x.distances.count(n)
        y.distances.count(n)
      }
      case _ => Unit
    }

    val metrics = annotated.map { ann =>
      val structure = worstStructure(ann.index)

      IndexMetric(
        index                     = new String(ann.index),
        source                    = if (ann.existing) "existing" else "novel",
        min_mismatches            = ann.distances.headOption.map { case (distance, count) => distance }.getOrElse(ann.index.length),
        indices_at_min_mismatches = ann.distances.headOption.map { case (distance, count) => count.toInt }.getOrElse(0),
        gc                        = SequenceUtil.calculateGc(ann.index),
        longest_homopolymer       = Sequences.longestHomopolymer(new String(ann.index))._2,
        worst_structure_seq       = structure.map(_.sequence()),
        worst_structure_dbn           = structure.map(_.structure()),
        worst_structure_delta_g   = structure.map(_.deltaG)
      )
    }

    Metric.write(output, metrics)
  }
}
