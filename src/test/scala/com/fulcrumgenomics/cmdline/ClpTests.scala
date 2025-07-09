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

package com.fulcrumgenomics.cmdline

import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Metric

import java.nio.file.Path

/* This a silly CLP. */
@clp(group=ClpGroups.Utilities, description="A test class")
class TestClp
(
  @arg(flag='e', doc="If set, exit with this code.")    val exitCode: Option[Int],
  @arg(flag='m', doc="If set, fail with this message.") val message: Option[String],
  @arg(flag='p', doc="Print this message.")             val printMe: Option[String],
  @arg(flag='i', doc="Write the tool information.")     val infoPath: Option[Path] = None
) extends FgBioTool {
  override def execute(): Unit = {
    (exitCode, message) match {
      case (Some(ex), Some(msg)) => fail(ex, msg)
      case (Some(ex), None     ) => fail(ex)
      case (None,     Some(msg)) => fail(msg)
      case (None,     None     ) => printMe.foreach(println)
    }
    this.infoPath.foreach { path =>
      this.toolInfo.foreach { info => Metric.write(path, info) }
    }
  }
}

/** Some basic test for the CLP classes. */
class ClpTests extends UnitSpec {
  "FgBioMain" should "find a CLP and successfully set it up and execute it" in {
    new FgBioMain().makeItSo("TestClp --print-me=hello".split(' ')) shouldBe 0
  }

  it should "fail with the provided exit code" in {
    new FgBioMain().makeItSo("TestClp -e 7".split(' ')) shouldBe 7
    new FgBioMain().makeItSo("TestClp --exit-code=5".split(' ')) shouldBe 5
    new FgBioMain().makeItSo("TestClp --exit-code=9 --message=FailBabyFail".split(' ')) shouldBe 9
    new FgBioMain().makeItSo("TestClp --message=FailBabyFail".split(' ')) should not be 0
  }

  it should "fail and print usage" in {
    new FgBioMain().makeItSo("SomeProgram --with-args=that-dont-exist".split(' ')) should not be 0
  }

  it should "provide command line information to the tool" in {
    val tmpPath = makeTempFile("tool_info.", ".tab")
    new FgBioMain().makeItSo(args=s"--async-io true TestClp -i $tmpPath".split(' ')) shouldBe 0
    val metrics = Metric.read[FgBioToolInfo](tmpPath)
    metrics should have length 1
    val metric = metrics.head

    // args/commandLineWithoutDefaults
    metric.commandLineWithDefaults should include ("fgbio --async-io true") // argument set _prior_ to the tool name
    metric.commandLineWithDefaults should not include ("fgbio --compression 5") // default argument _prior_ to the tool name
    metric.commandLineWithoutDefaults should include ("TestClp")
    metric.commandLineWithoutDefaults should include (s"-i $tmpPath")
    metric.commandLineWithoutDefaults should not include ("--print-me :none:") // a default argument
    metric.commandLineWithoutDefaults shouldBe s"fgbio --async-io true TestClp -i $tmpPath"

    // commandLineWithDefaults
    metric.commandLineWithDefaults should include ("fgbio --async-io true") // argument set _prior_ to the tool name
    metric.commandLineWithDefaults should include (" --compression 5") // default argument _prior_ to the tool namee
    metric.commandLineWithDefaults should include ("TestClp")
    metric.commandLineWithDefaults should include (s"--info-path $tmpPath")
    metric.commandLineWithDefaults should include ("--print-me :none:") // a default argument
    metric.commandLineWithDefaults shouldBe f"fgbio --async-io true --compression 5 --tmp-dir ${FgBioCommonArgs.args.tmpDir} --log-level Info --sam-validation-stringency SILENT TestClp --info-path ${tmpPath} --exit-code :none: --message :none: --print-me :none:"

    // description
    metric.description shouldBe "A test class"

    // version
    // Since the JAR has not been backaged yet, we cannot get the implementation version.
    metric.version shouldBe "null" // not set in scalatest, so unable to test
  }
}
