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

  "CopyUmiFromReadName" should "copy the UMI from a read name" in {
    val builder = new SamBuilder()

    builder.addFrag(name="1:AAAA", unmapped=true)
    builder.addFrag(name="1:2:CCCC", unmapped=true)
    builder.addFrag(name="1:2:3:GGGG", unmapped=true)
    builder.addFrag(name="blah:AAAA+CCCC", unmapped=true)

    val out  = makeTempFile("test.", ".bam")
    val tool = new CopyUmiFromReadName(input=builder.toTempFile(), output=out)
    executeFgbioTool(tool)

    val umis = readBamRecs(out).map(rec => rec[String]("RX"))
    umis should contain theSameElementsInOrderAs Seq("AAAA", "CCCC", "GGGG", "AAAA-CCCC")
  }
}
