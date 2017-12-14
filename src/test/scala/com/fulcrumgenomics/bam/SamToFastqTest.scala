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

package com.fulcrumgenomics.bam

import java.nio.file.Files

import com.fulcrumgenomics.FgBioDef.{DirPath, FilePath}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.fastq.{FastqRecord, FastqSource}
import com.fulcrumgenomics.testing.SamBuilder._
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import org.scalatest.OptionValues

class SamToFastqTest extends UnitSpec with OptionValues {

  private def addTags(rec: SamRecord, b5: Option[String] = None, q5: Option[String] = None, b3: Option[String] = None, q3: Option[String] = None): Unit = {
    b5.foreach { value => rec(SamToFastq.FivePrimeBasesTag)      = value}
    q5.foreach { value => rec(SamToFastq.FivePrimeQualitiesTag)  = value}
    b3.foreach { value => rec(SamToFastq.ThreePrimeBasesTag)     = value}
    q3.foreach { value => rec(SamToFastq.ThreePrimeQualitiesTag) = value}
  }

  private object FastqResults {
    def apply(path: DirPath): FastqResults = {
      FastqResults(
        frags = FastqSource(PathUtil.pathTo(path + ".fragments.fq.gz")).toSeq,
        r1s   = FastqSource(PathUtil.pathTo(path + ".R1.fq.gz")).toSeq,
        r2s   = FastqSource(PathUtil.pathTo(path + ".R2.fq.gz")).toSeq
      )
    }
  }

  private case class FastqResults(frags: Seq[FastqRecord], r1s: Seq[FastqRecord], r2s: Seq[FastqRecord])

  private def checkFastqRecord(rec: FastqRecord, name: String, bases: String, quals: String): Unit = {
    rec.name shouldBe name
    rec.bases shouldBe bases
    rec.quals shouldBe quals
  }

  private def toOutput: FilePath = {
    val dir = Files.createTempDirectory("SamToFastq")
    dir.toFile.deleteOnExit()
    PathUtil.pathTo(dir.toAbsolutePath.toString, "fastqs")
  }

  "SamToFastq" should "convert fragment records" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname), readLength=10)
    builder.addFrag(name="q1", bases="A"*10, quals="I"*10, unmapped=true).foreach { rec => addTags(rec=rec, b5=Some("A"), q5=Some("5"), b3=Some("C"), q3=Some("3")) }
    builder.addFrag(name="q2", bases="A"*10, quals="I"*10, start=10, strand=Plus).foreach { rec => addTags(rec=rec, b5=Some("C"), q5=Some("6"), b3=Some("G"), q3=Some("4")) }
    builder.addFrag(name="q3", bases="A"*10, quals="I"*5 + "J"*5, start=10, strand=Minus).foreach { rec => addTags(rec=rec, b5=Some("G"), q5=Some("7"), b3=Some("T"), q3=Some("5")) }

    val output = toOutput
    new SamToFastq(input=builder.toTempFile(), output=output).execute()
    
    val results = FastqResults(output)
    results.frags should have size 3
    results.r1s shouldBe 'empty
    results.r2s shouldBe 'empty

    checkFastqRecord(rec=results.frags(0), name="q1", bases="A"*11 + "C",       quals="5" + "I"*10 + "3")
    checkFastqRecord(rec=results.frags(1), name="q2", bases="C" + "A"*10 + "G", quals="6" + "I"*10 + "4")
    checkFastqRecord(rec=results.frags(2), name="q3", bases="G" + "T"*10 + "T", quals="7" + "J"*5 + "I"*5 + "5")
  }

  it should "convert paired records" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname), readLength=10)
    // unmapped
    builder.addPair(name="q1", bases1="A"*10, quals1="I"*10, unmapped1=true, bases2="A"*10, quals2="I"*10, unmapped2=true).foreach { rec =>
      addTags(rec=rec, b5=Some("A"), q5=Some("5"), b3=Some("C"), q3=Some("3"))
    }
    // FR
    builder.addPair(name="q2", bases1="A"*10, quals1="I"*5 + "J"*5, start1=10, bases2="A"*10, quals2="I"*5 + "J"*5, start2=10).foreach { rec =>
      addTags(rec=rec, b5=Some("C"), q5=Some("6"), b3=Some("G"), q3=Some("4")) }
    // RF
    builder.addPair(name="q3", bases1="A"*10, quals1="I"*5 + "J"*5, start1=10, strand1=Minus, bases2="A"*10, quals2="I"*5 + "J"*5, start2=10, strand2=Plus).foreach { rec =>
      addTags(rec=rec, b5=Some("G"), q5=Some("7"), b3=Some("T"), q3=Some("5")) }

    val output = toOutput
    new SamToFastq(input=builder.toTempFile(), output=output).execute()

    val results = FastqResults(output)
    results.frags shouldBe 'empty
    results.r1s should have size 3
    results.r2s should have size 3

    checkFastqRecord(rec=results.r1s(0), name="q1", bases="A"*11 + "C",       quals="5" + "I"*10 + "3")
    checkFastqRecord(rec=results.r1s(1), name="q2", bases="C" + "A"*10 + "G", quals="6" + "I"*5 + "J"*5 + "4") // forward
    checkFastqRecord(rec=results.r1s(2), name="q3", bases="G" + "T"*10 + "T", quals="7" + "J"*5 + "I"*5 + "5") // reverse

    checkFastqRecord(rec=results.r2s(0), name="q1", bases="A" + "T"*10 + "C", quals="5" + "I"*10 + "3")
    checkFastqRecord(rec=results.r2s(1), name="q2", bases="C" + "T"*10 + "G", quals="6" + "J"*5 + "I"*5 + "4") // reverse
    checkFastqRecord(rec=results.r2s(2), name="q3", bases="G" + "A"*10 + "T", quals="7" + "I"*5 + "J"*5 + "5") // forward
  }

  it should "not prepend or append if no tags are given" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Queryname), readLength=10)
    builder.addFrag(name="q1", bases="A"*10, quals="I"*10, unmapped=true)
    builder.addPair(name="q2", bases1="A"*10, quals1="I"*10, unmapped1=true, bases2="A"*10, quals2="I"*10, unmapped2=true)

    val output = toOutput
    new SamToFastq(input=builder.toTempFile(), output=output).execute()

    val results = FastqResults(output)
    results.frags should have size 1
    results.r1s should have size 1
    results.r2s should have size 1

    checkFastqRecord(rec=results.frags(0), name="q1", bases="A"*10, quals="I"*10)
    checkFastqRecord(rec=results.r1s(0),   name="q2", bases="A"*10, quals="I"*10)
    checkFastqRecord(rec=results.r2s(0),   name="q2", bases="T"*10, quals="I"*10)
  }

  it should "sort the input to queryname" in {
    val builder = new SamBuilder(sort=Some(SamOrder.Coordinate), readLength=10)
    builder.addFrag(name="q1", bases="A"*10, quals="I"*10, unmapped=true)
    builder.addFrag(name="q2", bases="A"*10, quals="I"*10, start=100)
    builder.addFrag(name="q3", bases="A"*10, quals="I"*10, start=50)

    val output = toOutput
    new SamToFastq(input=builder.toTempFile(), output=output).execute()

    builder.map(_.name).toSeq should contain theSameElementsInOrderAs Seq("q3", "q2", "q1")

    val results = FastqResults(output)
    results.frags should have size 3
    results.r1s shouldBe 'empty
    results.r2s shouldBe 'empty

    checkFastqRecord(rec=results.frags(0), name="q1", bases="A"*10, quals="I"*10)
    checkFastqRecord(rec=results.frags(1), name="q2", bases="A"*10, quals="I"*10)
    checkFastqRecord(rec=results.frags(2), name="q3", bases="A"*10, quals="I"*10)
  }
}
