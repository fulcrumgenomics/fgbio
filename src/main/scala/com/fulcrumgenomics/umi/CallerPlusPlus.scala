/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics
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

package com.fulcrumgenomics.umi

import java.io.{BufferedReader, BufferedWriter, InputStreamReader, OutputStreamWriter}
import java.nio.file.Paths

import com.fulcrumgenomics.FgBioDef.{FgBioEnum, FilePath, forloop}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.sopt.util._
import com.fulcrumgenomics.umi.ConsensusCaller.Base
import com.fulcrumgenomics.umi.UmiConsensusCaller.SourceRead
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.util.NumericTypes.PhredScore
import enumeratum.EnumEntry
import htsjdk.samtools.SAMUtils

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable.ListBuffer


object CallerPlusPlus {

  case class Query(name: String, sequences: Seq[String]) {
    def fastaString: String = ">" + name + "\n" + sequences.mkString("\n") + "\n"
  }

  case class Result(name: String, consensus: String, msa: Seq[String] = Seq.empty)
}


sealed trait AlignmentAlgorithm extends EnumEntry {
  def value: Int
}
object AlignmentAlgorithm extends FgBioEnum[AlignmentAlgorithm] {
  def values: IndexedSeq[AlignmentAlgorithm] = findValues
  case object Local extends AlignmentAlgorithm { val value: Int = 0 }
  case object Global extends AlignmentAlgorithm { val value: Int = 1 }
  case object SemiGlobal extends AlignmentAlgorithm { val value: Int = 2 }
}


class CallerPlusPlus(val callerpp: FilePath = Paths.get("callerpp"),
                     val matchScore: Int = 5,
                     val mismatchScore: Int = -4,
                     val gapScore: Int = -8,
                     val algorithm: AlignmentAlgorithm = AlignmentAlgorithm.Local,
                     val pairwiseAlignment: Boolean = false) extends LazyLogging with AutoCloseable {

  import CallerPlusPlus._

  private val args: Seq[String] = {
    val buffer = ListBuffer[Any]()
    buffer.append(callerpp)
    buffer.append("-A", matchScore)
    buffer.append("-B", mismatchScore)
    buffer.append("-O", gapScore)
    buffer.append("-a", algorithm.value)
    buffer.append("-m") // output the MSA
    buffer.map(_.toString)
  }

  //logger.info("Proccess args: " + args.mkString(" "))

  private val process       = new ProcessBuilder(args:_*).start()
  private val errorPipe     = Io.pipeStream(process.getErrorStream, s => logger.debug("callerpp: ", s))
  private val toCallerpp    = new BufferedWriter(new OutputStreamWriter(process.getOutputStream))
  private val fromCallerpp  = new BufferedReader(new InputStreamReader(process.getInputStream))

  def map(query: Query): Result = map(Seq(query)).head

  def map(queries: Seq[Query]): Seq[Result] = if (queries.isEmpty) Seq.empty else {

    // write them to the input of callerpp
    queries.foreach { query => toCallerpp.write(query.fastaString) }
    signalCallerpp()

    // read the results
    queries.map { query =>
      // Read the name
      val name = this.fromCallerpp.readLine() match {
        case null  => throw new IllegalStateException(s"callerpp failed on:\n${query.sequences.mkString("\n")}")
        case _name =>
          assert(_name.startsWith(">"), s"Invalid name: ${_name}")
          _name.drop(1)
      }
      assert(name == query.name, s"Query and Result are out of order: Query=${query.name} Result=$name")

      // Read the consensus
      val consensus = this.fromCallerpp.readLine()

      // Read the MSA (NB: an extra for the consensus sequence)
      val msa = {
        // read in the alignment strings for the MSA
        val alignmentStrings = Range.inclusive(0, query.sequences.length).map { _ => this.fromCallerpp.readLine() }

        // convert the alignment strings to bases
        val alignmentBases = alignmentStrings.map(_.replaceAll("[ -.]", "")).toArray // just keep the bases

        // NB: the consensus sequence is always at the end
        val sequences = query.sequences :+ consensus

        // verify that the input sequences match the MSA in order
        assert(alignmentBases.zip(sequences).forall { case (a, b) => a == b }, "Bug: msa out of order")

        // move the consensus to the first sequence
        alignmentStrings.last +: alignmentStrings.dropRight(1)
      }
      assert(msa.length == query.sequences.length + 1)

      // Create the result
      Result(name=name, consensus=consensus, msa=msa)
    }.toList
  }

  private def signalCallerpp(): Unit = {
    this.toCallerpp.flush()
    this.toCallerpp.newLine()
    this.toCallerpp.flush()
  }

  override def close(): Unit = {
    this.toCallerpp.close()
    this.fromCallerpp.close()
    this.process.destroy()
    this.errorPipe.close(2000)
  }
}


private case class ConsensusMsaInfo(consensusMsaBase: Char, consensusBase: Base, consensusQual: PhredScore, baseCounter: SimpleCounter[Char] = new SimpleCounter[Char]()) {

  private def toChar(count: Long): Char = if (count < 10) s"$count".head else '+'

  val depth: Char  = toChar(baseCounter.map { case (b, count) => if (b == '-') 0 else count }.sum)
  val errors: Char = toChar(baseCounter.map { case (b, count) => if (b == consensusMsaBase || b == '-') 0 else count }.sum)
  def countOf(base: Char): Char = toChar(this.baseCounter.countOf(base))
}


case object KBLDGRN extends TermCode {
  override val code: String = "\u001B[1m\u001B[32m"
}

class CallerPlusPlusConsensusCaller(val options: VanillaUmiConsensusCallerOptions = new VanillaUmiConsensusCallerOptions(),
                                    msaDebug: Boolean = false
                                   ) extends ConsensusCallerTrait with LazyLogging {
  import CallerPlusPlus._

  private val callerPlusPlus  = new CallerPlusPlus(pairwiseAlignment=true)
  private val consensusCaller = new ConsensusCaller(
    errorRatePreLabeling  = this.options.errorRatePreUmi,
    errorRatePostLabeling = this.options.errorRatePostUmi
  )
  private val NoCall: Byte = 'N'.toByte
  private val NotEnoughReadsQual: PhredScore = 0.toByte // Score output when masking to N due to insufficient input reads
  private val TooLowQualityQual: PhredScore = 2.toByte  // Score output when masking to N due to too low consensus quality
  private val MismatchQual: PhredScore = 3.toByte // Score output when the POA consensus call mismatches the vanilla consensus call

  def consensusCall(reads: Seq[SourceRead]): Option[VanillaConsensusRead] = if (reads.isEmpty) None else {
    // Create the query
    val query    = Query(name=reads.head.id, sequences=reads.map(_.baseString))

    // Run the partial order aligner and generate a consensus sequence
    val result   = callerPlusPlus.map(query)

    // get the most likely consensus bases and qualities
    val builder         = this.consensusCaller.builder()
    val consensusLength = result.consensus.length
    val consensusBases  = new Array[Base](consensusLength)
    val consensusQuals  = new Array[PhredScore](consensusLength)
    val consensusDepths = new Array[Short](consensusLength)
    val consensusErrors = new Array[Short](consensusLength)
    val consensusMsa    = result.msa.head
    val readsInMsa      = result.msa.drop(1)
    var consensusOffset = 0
    val readOffsets = Array.fill(n=readsInMsa.length)(0)
    forloop(from = 0, until = consensusMsa.length) { msaIndex =>
      if (consensusMsa(msaIndex) != '-') {
        forloop(from = 0, until = readsInMsa.length) { readIndex =>
          val readInMsa = readsInMsa(readIndex)
          if (readInMsa(msaIndex) != '-') {
            val read       = reads(readIndex)
            val readOffset = readOffsets(readIndex)
            val readBase   = read.bases(readOffset)
            val readQual   = read.quals(readOffset)
            builder.add(readBase, readQual)
            readOffsets(readIndex) += 1
          }
        }

        val cBase = result.consensus(consensusOffset).toByte
        val (rawBase, rawQual) = builder.call()
        val (base, qual) = {
          if (builder.contributions < this.options.minReads)       (NoCall, NotEnoughReadsQual)
          else if (rawQual < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual)
          else if (rawBase != cBase)                               (cBase,  MismatchQual)
          else (rawBase, rawQual)
        }

        consensusBases(consensusOffset) = base
        consensusQuals(consensusOffset) = qual

        // Generate the values for depth and count of errors
        val depth  = builder.contributions
        val errors = if (rawBase == NoCall) depth else depth - builder.observations(rawBase)
        consensusDepths(consensusOffset) = if (depth  > Short.MaxValue) Short.MaxValue else depth.toShort
        consensusErrors(consensusOffset) = if (errors > Short.MaxValue) Short.MaxValue else errors.toShort

        consensusOffset += 1
        builder.reset()
      }
    }

    assert(result.consensus.length == consensusBases.length)
    /*
    logger.info(s"name: ${result.name}")
    logger.info(s"consensus 1: '${result.consensus}'")
    logger.info(s"consensus 2: '${new String(consensusBases)}")
    logger.info(s"depth: " + consensusDepths.map(_.toString).mkString(","))
    logger.info(s"errors: " + consensusErrors.map(_.toString).mkString(","))
    result.msa.foreach { msaLine => logger.info("msa: " + msaLine) }
    */
    if (msaDebug) {
      println(KBLD(KGRN(s">${result.name}")))

      var consensusIndex = 0
      val trimmedMsas = result.msa.drop(1).map { msa =>
        msa.reverseIterator.dropWhile(_ == '-').mkString.reverse
      }
      val infos = consensusMsa.indices.map { idx =>
        consensusMsa(idx) match {
          case '-' =>
            ConsensusMsaInfo(consensusMsaBase='-', consensusBase='-'.toByte, consensusQual='-'.toByte)
          case base =>
            val consensusBase  = consensusBases(consensusIndex)
            val consensusQual  = consensusQuals(consensusIndex)
            val baseCounter    = SimpleCounter[Char](trimmedMsas.flatMap(m => if (idx < m.length) Some(m(idx)) else None))
            consensusIndex += 1
            ConsensusMsaInfo(consensusMsaBase=base, consensusBase=consensusBase, consensusQual=consensusQual, baseCounter=baseCounter)
        }
      }
      val defaultPadding = 30
      def printlnPad(str: String, padding: Int): Unit = println((" " * padding) + str)
      // called consensus base
      print(KBLDGRN("[bases]"))
      printlnPad(KBLDGRN(infos.map(_.consensusBase.toChar).mkString), padding=defaultPadding-7)
      // called consensus quality
      print(KBLDGRN("[quals]"))
      printlnPad(KBLDGRN(infos.map(info => SAMUtils.phredToFastq(info.consensusQual.toInt)).mkString), padding=defaultPadding-7)
      // counts of ACGTN-s
      "ACGTN-".foreach { base: Char =>
        val colorMap = infos.map(info => if (info.consensusBase == base) KBLDGRN else KBLDRED)
        val label    = s"[$base]"
        val bases    = infos.zip(colorMap).map { case (info, colorFn) =>
          info.countOf(base).toString match {
            case "0"   => KGRN("-")
            case count => colorFn(count)
          }
        }.mkString
        print(KGRN(label))
        printlnPad(bases, padding=defaultPadding-label.length)
      }
      def colorDepth(count: Char): String = if (count == '0') KBLDRED("0") else KBLDGRN(count.toString)
      def colorError(count: Char): String = if (count == '0') KBLDGRN("0") else KBLDRED(count.toString)
      // depths ('+' means > 9)
      print(KBLDGRN("[depth]"))
      printlnPad(infos.map(_.depth).map(colorDepth).mkString, padding=defaultPadding-7)
      // errors ('+' means > 9)
      print(KBLDGRN("[errors]"))
      printlnPad(infos.map(_.errors).map(colorError).mkString, padding=defaultPadding-8)

      // MSA
      val sequences: Seq[String] = result.consensus +: query.sequences
      forloop(from=0, until=sequences.length) { idx =>
        val name      = KCYN(f"type=${if (0 == idx) "C" else s"$idx"}%-3s")
        val inLength  = sequences(idx).length
        val msa       = colorMsa(result.msa(idx), consensusMsa)
        val outLength = result.msa(idx).count(_ != '-')
        val lengthCode: String => String = str => if (inLength == outLength) KBLDGRN(str) else KBLDRED(str)
        val inLengthStr  = lengthCode(f"in=$inLength%3d")
        val outLengthStr = lengthCode(f"out=$outLength%3d")
        val remaining    = KBLDRED(sequences(idx).substring(outLength))
        println(f"[$inLengthStr] [$outLengthStr] [$name] $msa$remaining")
      }
      //print("Waiting on <enter>...")
      //scala.io.StdIn.readLine()
    }

    val read = VanillaConsensusRead(
      id     = result.name,
      bases  = consensusBases,
      quals  = consensusQuals,
      depths = consensusDepths,
      errors = consensusErrors
    )
    Some(read)
  }

  private def colorMsa(msa: String, consensusMsa: String): String = {
    var startIdx = 0
    val buffer = new ListBuffer[(Int, String)]()

    // 0 - match
    // 1 - mismatch
    // 2 - dash
    var state = 0
    forloop(from=0, until=msa.length) { endIndex =>
      (msa(endIndex), consensusMsa(endIndex)) match {
        case ('-', _) =>
          if (state != 2) {
            buffer += ((state, msa.substring(startIdx, endIndex)))
            state = 2
            startIdx = endIndex
          }
        case (m, c) if m != c =>
          if (state != 1) {
            buffer += ((state, msa.substring(startIdx, endIndex)))
            state = 1
            startIdx = endIndex
          }
        case (_, _) =>
          if (state != 0) {
            buffer += ((state, msa.substring(startIdx, endIndex)))
            state = 0
            startIdx = endIndex
          }
      }
    }
    buffer += ((state, msa.substring(startIdx)))

    assert(msa.length == buffer.map(_._2.length).sum, s"buffer: ${buffer.toList}")
    buffer.map { case (st, subMsa) =>
      st match {
        case 0     => KBLDGRN(subMsa)
        case 1 | 2 => KBLDRED(subMsa)
      }
    }.mkString
  }
}
