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

package com.fulcrumgenomics.bam

import java.nio.file.Paths

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import org.apache.commons.math3.stat.regression.SimpleRegression

import scala.util.Random

/**
  * Tests for RandomizeBam
  */
class RandomizeBamTest extends UnitSpec {
  private val bam = Paths.get("src/test/resources/com/fulcrumgenomics/bam/200reads.bam")

  /** Slurps read names with /1 or /2 suffixes, in order, from a BAM file. */
  def slurp(p: PathToBam): IndexedSeq[String] = readBamRecs(p).map(_.id)

  "RandomizeBam" should "truly randomize the order of reads in a file in non-query-grouped mode" in {
    val out1 = makeTempFile("random1.", ".bam")
    val out2 = makeTempFile("random1.", ".bam")
    val random = new Random(1)
    Seq(out1, out2).foreach(out => new RandomizeBam(input=bam, output=out, seed=random.nextInt(), queryGroup=false).execute())

    val in = slurp(bam)
    val o1 = slurp(out1)
    val o2 = slurp(out2)

    val inPos = in.zipWithIndex.map(_._2)
    val o1Pos = o1.map(name => in.indexOf(name))
    val o2Pos = o2.map(name => in.indexOf(name))

    Seq((inPos, o1Pos), (inPos, o2Pos), (o1Pos, o2Pos)).foreach { case (is1, is2) => {
      val regression = new SimpleRegression()
      is1.zip(is2).foreach(pair => regression.addData(pair._1, pair._2))
      regression.regress().getRSquared should be < 0.05
    }}
  }

  it should "randomize the order of reads in a file in query-grouped mode" in {
    val out1 = makeTempFile("random1.", ".bam")
    val out2 = makeTempFile("random1.", ".bam")
    val random = new Random(7)
    Seq(out1, out2).foreach(out => new RandomizeBam(input=bam, output=out, seed=random.nextInt(), queryGroup=true).execute())

    val in = slurp(bam)
    val o1 = slurp(out1)
    val o2 = slurp(out2)

    val inPos = in.zipWithIndex.map(_._2)
    val o1Pos = o1.map(name => in.indexOf(name))
    val o2Pos = o2.map(name => in.indexOf(name))

    Seq((inPos, o1Pos), (inPos, o2Pos), (o1Pos, o2Pos)).foreach { case (is1, is2) => {
      val regression = new SimpleRegression()
      is1.zip(is2).foreach(pair => regression.addData(pair._1, pair._2))
      regression.regress().getRSquared should be < 0.05
    }}

    // Additionally validate that the outputs are query-grouped
    Seq(o1, o2).foreach(names => {
      names.grouped(2).foreach { case Seq(r1, r2) => {
        r1.substring(0, r1.length-2) shouldEqual r2.substring(0, r2.length-2)
      }}
    })
  }
}
