/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics
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

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import org.scalatest.OptionValues

class CopyUmiFromReadNameTest extends UnitSpec with OptionValues {

  private case class Result(name: String, umi: String)

  /** Runs CopyUmiFromReadName using the given read names returning the output read names and UMIs. */
  private def run(names: Iterable[String], removeUmi: Boolean, umiDelimiter: Option[Char]): IndexedSeq[Result] = {
    // build the reads
    val builder = new SamBuilder()
    names.foreach { name => builder.addFrag(name=name, unmapped=true) }

    // run the tool
    val out  = makeTempFile("test.", ".bam")
    val tool = new CopyUmiFromReadName(input=builder.toTempFile(), output=out, removeUmi=removeUmi, umiDelimiter=umiDelimiter)
    executeFgbioTool(tool)

    // slurp the results
    val recs = readBamRecs(out)
    recs.length shouldBe builder.size
    recs.map { rec => Result(name=rec.name, umi=rec[String](ConsensusTags.UmiBases)) }
  }

  "CopyUmiFromReadName" should "copy the UMI from a read name" in {
    val names   = Seq("1:AAAA", "1:2:CCCC", "1:2:3:GGGG", "blah:AAAA-CCCC")
    val results = run(names=names, removeUmi=false, umiDelimiter=None)
    results.map(_.name) should contain theSameElementsInOrderAs names
    results.map(_.umi) should contain theSameElementsInOrderAs Seq("AAAA", "CCCC", "GGGG", "AAAA-CCCC")
  }

  it should "remove the UMI from the read name when --remove-umi=true" in {
    val names   = Seq("1:AAAA", "1:2:CCCC", "1:2:3:GGGG", "blah:AAAA-CCCC")
    val results = run(names=names, removeUmi=true, umiDelimiter=None)
    results.map(_.name) should contain theSameElementsInOrderAs Seq("1", "1:2", "1:2:3", "blah")
    results.map(_.umi) should contain theSameElementsInOrderAs Seq("AAAA", "CCCC", "GGGG", "AAAA-CCCC")
  }

  it should "update the UMI delimiter in the read name when --umi-delimiter=+" in {
    val names   = Seq("1:AAAA", "1:2:CCCC", "1:2:3:GGGG", "blah:AAAA+CCCC")
    val results = run(names=names, removeUmi=true, umiDelimiter=Some('+'))
    results.map(_.name) should contain theSameElementsInOrderAs Seq("1", "1:2", "1:2:3", "blah")
    results.map(_.umi) should contain theSameElementsInOrderAs Seq("AAAA", "CCCC", "GGGG", "AAAA-CCCC")
  }
}
