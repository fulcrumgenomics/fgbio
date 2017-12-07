/*
 * The MIT License
 *
 * Copyright (c) 2017 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.commons.util.SimpleCounter
import com.fulcrumgenomics.umi.ConsensusCaller.Base
import com.fulcrumgenomics.umi.UmiConsensusCaller.SourceRead
import com.fulcrumgenomics.util.NumericTypes.PhredScore

import scala.collection.mutable.ArrayBuffer

class GraphBaseConsensusCaller(msaCommand: String = DuplexConsensusCaller.MsaCommand,
                               fallbackToVanillaCaller: Boolean = false,
                               options: VanillaUmiConsensusCallerOptions = new VanillaUmiConsensusCallerOptions())
  extends MultipleSequenceAlignmentConsensusCaller(msaCommand=msaCommand, fallbackToVanillaCaller=fallbackToVanillaCaller, options=options) {


  override protected def call(aligned: Seq[SourceRead]): VanillaConsensusRead = {
    val expectedLength = aligned.head.length
    require(aligned.forall(_.length == expectedLength))

    // One extra node for the source, and one extra node for the sink
    val nodes = Array.range(0, expectedLength+2).map { offset => Node(offset=offset-1)}

    // add the edges
    aligned.map(_.bases).foreach { bases =>
      // NB: i an j are 1-based relative to the bases, and 0-based relative to the nodes
      var i = 0
      while (i < nodes.length - 1) { // all but the last node (sink)
        var j = i + 1
        while (j-1 < bases.length && bases(j-1) == Gap) {
          j = j + 1
        }
        Node.connect(nodes(i), nodes(j))
        i = j
      }
    }

    // sanity checks
    require(nodes.head.incoming.isEmpty)
    require(nodes.head.outgoing.map(_._2).sum == aligned.length)
    require(nodes.last.outgoing.isEmpty)
    require(nodes.last.incoming.map(_._2).sum == aligned.length)
    nodes.foreach { node =>
      require(node.incoming.map(_._2).sum <= aligned.length)
      require(node.outgoing.map(_._2).sum <= aligned.length)
      node.incoming.foreach{ case (inc, count) =>
        require(inc.outgoingCount(node) == count)
      }
      node.outgoing.foreach { case (out, count) =>
        require(out.incomingCount(node) == count)
      }
    }

    // Get the consensus alignment, returned as a list of offsets
    val offsets = {
      import scala.collection.mutable
      val from   = new mutable.HashMap[Node, Node]() // cur -> predecessor
      val scores = new mutable.HashMap[Node, Long]()

      scores(nodes.head) = 0

      nodes.drop(1).foreach { node =>
        val incoming = node.incoming.map { case (prevNode, count) =>
          // the score is the sum of the previous node's score and the count of the edge
          (prevNode, count + scores(prevNode))
        }
        val (bestIncoming, bestScore) = incoming.maxBy(_._2)
        from(node)   = bestIncoming
        scores(node) = bestScore
      }

      // traceback
      var curNode = nodes.last
      val bestPath = new ArrayBuffer[Node]()
      bestPath += curNode
      while (from.contains(curNode)) {
        val prevNode = from(curNode)
        bestPath += prevNode
        curNode = prevNode
      }
      bestPath.reverse
        .map(_.offset)
        .filter { offset => 0 <= offset && offset < expectedLength}
    }

    // Call the consensus read
    {
      val consensusLength = offsets.length
      val consensusBases  = new Array[Base](consensusLength)
      val consensusQuals  = new Array[PhredScore](consensusLength)
      val consensusDepths = new Array[Short](consensusLength)
      val consensusErrors = new Array[Short](consensusLength)

      val builder = this.caller.builder()
      offsets.zipWithIndex.foreach { case (offset, consensusOffset) =>
        // Add the evidence
        aligned.foreach { read =>
          val base = read.bases(offset)
          val qual = read.quals(offset)
          if (base != NoCall && base != Gap)  builder.add(base=base, qual=qual)
        }

        // Call the consensus and do any additional filtering
        val (rawBase, rawQual) = builder.call()
        val (base, qual) = {
          if (builder.contributions < this.options.minReads)       (NoCall, NotEnoughReadsQual)
          else if (rawQual < this.options.minConsensusBaseQuality) (NoCall, TooLowQualityQual)
          else (rawBase, rawQual)
        }

        consensusBases(consensusOffset) = base
        consensusQuals(consensusOffset) = qual

        // Generate the values for depth and count of errors
        val depth  = builder.contributions
        val errors = if (rawBase == NoCall) depth else depth - builder.observations(rawBase)
        consensusDepths(consensusOffset) = if (depth  > Short.MaxValue) Short.MaxValue else depth.toShort
        consensusErrors(consensusOffset) = if (errors > Short.MaxValue) Short.MaxValue else errors.toShort

        // Get ready for the next pass
        builder.reset()
      }
      VanillaConsensusRead(
        id     = aligned.head.id,
        bases  = consensusBases,
        quals  = consensusQuals,
        depths = consensusDepths,
        errors = consensusErrors
      )
    }
  }
}

private object Node {
  def connect(from: Node, to:Node): Unit = {
    from.addSuccessor(to)
    to.addPredecessor(from)
  }
}

private case class Node(offset: Int) {
  private val _incoming = new SimpleCounter[Node]
  private val _outgoing = new SimpleCounter[Node]

  def incoming: Seq[(Node, Long)] = this._incoming.toSeq
  def outgoing: Seq[(Node, Long)] = this._outgoing.toSeq

  def incomingCount(node: Node): Long = this._incoming.countOf(node)
  def outgoingCount(node: Node): Long = this._outgoing.countOf(node)

  def addSuccessor(successor: Node): Unit = {
    this._outgoing.count(successor)
  }

  def addPredecessor(predecessor: Node): Unit = {
    this._incoming.count(predecessor)
  }
}