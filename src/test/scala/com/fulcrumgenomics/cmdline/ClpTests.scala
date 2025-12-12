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

import com.fulcrumgenomics.FgBioDef.DirPath
import com.fulcrumgenomics.bam.api.{SamSource, SamWriter}
import com.fulcrumgenomics.commons.util.{LogLevel, Logger}
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.{Io, Metric}
import htsjdk.samtools.ValidationStringency

import java.nio.file.{Path, Paths}

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

  /** Captures the state of all global variables that FgBioCommonArgs mutates. */
  private case class GlobalStateSnapshot(
    samSourceAsyncIo: Boolean,
    samWriterAsyncIo: Boolean,
    samWriterCompression: Int,
    ioCompression: Int,
    ioTmpDir: DirPath,
    loggerLevel: LogLevel,
    samSourceValidation: ValidationStringency,
    javaIoTmpDir: String,
    fgBioArgs: FgBioCommonArgs
  )

  /** Captures the current state of all global variables. */
  private def captureGlobalState(): GlobalStateSnapshot = GlobalStateSnapshot(
    samSourceAsyncIo = SamSource.DefaultUseAsyncIo,
    samWriterAsyncIo = SamWriter.DefaultUseAsyncIo,
    samWriterCompression = SamWriter.DefaultCompressionLevel,
    ioCompression = Io.compressionLevel,
    ioTmpDir = Io.tmpDir,
    loggerLevel = Logger.level,
    samSourceValidation = SamSource.DefaultValidationStringency,
    javaIoTmpDir = System.getProperty("java.io.tmpdir"),
    fgBioArgs = FgBioCommonArgs.args
  )

  /** Restores the global state from a snapshot. */
  private def restoreGlobalState(snapshot: GlobalStateSnapshot): Unit = {
    SamSource.DefaultUseAsyncIo = snapshot.samSourceAsyncIo
    SamWriter.DefaultUseAsyncIo = snapshot.samWriterAsyncIo
    SamWriter.DefaultCompressionLevel = snapshot.samWriterCompression
    Io.compressionLevel = snapshot.ioCompression
    Io.tmpDir = snapshot.ioTmpDir
    Logger.level = snapshot.loggerLevel
    SamSource.DefaultValidationStringency = snapshot.samSourceValidation
    System.setProperty("java.io.tmpdir", snapshot.javaIoTmpDir)
    FgBioCommonArgs.args = snapshot.fgBioArgs
  }

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
    metric.commandLineWithDefaults shouldBe s"fgbio --async-io true --compression 5 --tmp-dir ${FgBioCommonArgs.args.tmpDir} --log-level Info --sam-validation-stringency SILENT --cram-ref-fasta :none: TestClp --info-path ${tmpPath} --exit-code :none: --message :none: --print-me :none:"

    // description
    metric.description shouldBe "A test class"

    // version
    // Since the JAR has not been backaged yet, we cannot get the implementation version.
    metric.version shouldBe "null" // not set in scalatest, so unable to test
  }

  "FgBioCommonArgs" should "store and retrieve cramRefFasta value" in {
    val refPath = Some(Paths.get("/path/to/reference.fasta"))
    val args = new FgBioCommonArgs(cramRefFasta = refPath)
    args.cramRefFasta shouldBe refPath
  }

  it should "default cramRefFasta to None" in {
    val args = new FgBioCommonArgs()
    args.cramRefFasta shouldBe None
  }

  it should "provide cramRefFasta as default to SamSource and SamWriter" in {
    val snapshot = captureGlobalState()
    try {
      // Set common args with None (testing with BAM, not CRAM)
      FgBioCommonArgs.args = new FgBioCommonArgs(cramRefFasta = None)

      // Create test data
      val builder = new SamBuilder()
      builder.addPair(name="q1", start1=100, start2=300)
      val bam = makeTempFile("test.", ".bam")

      // Write and read should use common arg default
      val writer = SamWriter(bam, builder.header)
      writer ++= builder.iterator
      writer.close()

      val reader = SamSource(bam)
      val records = reader.toSeq
      records.size shouldBe 2
      reader.close()
    } finally {
      restoreGlobalState(snapshot)
    }
  }

  it should "use cramRefFasta from common args when writing and reading CRAM files" in {
    val snapshot = captureGlobalState()
    try {
      // Create a small reference sequence (500bp chr1)
      val refSequence = "N" * 500  // Simple 500bp sequence of N's
      val refLen = refSequence.length

      // Create synthetic reference FASTA file
      val ref = makeTempFile("reference.", ".fasta")
      val refWriter = Io.toWriter(ref)
      refWriter.write(s">chr1\n")
      refWriter.write(refSequence + "\n")
      refWriter.close()

      // Create .fai index file for the reference
      val fai = Paths.get(ref.toString + ".fai")
      val faiWriter = Io.toWriter(fai)
      faiWriter.write(s"chr1\t$refLen\t6\t$refLen\t${refLen + 1}\n")  // name, length, offset, basesPerLine, bytesPerLine
      faiWriter.close()

      // Set common args with the reference path
      FgBioCommonArgs.args = new FgBioCommonArgs(cramRefFasta = Some(ref))

      // Create test data with custom sequence dictionary matching our small reference
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=refLen))
      val builder = new SamBuilder(sd = Some(dict))
      builder.addPair(name="q1", start1=100, start2=200)  // Coordinates within 500bp
      val cram = makeTempFile("test.", ".cram")

      // Write to CRAM without explicitly passing ref - should use cramRefFasta from common args
      val writer = SamWriter(cram, builder.header)
      writer ++= builder.iterator
      writer.close()

      // Read from CRAM without explicitly passing ref - should use cramRefFasta from common args
      val reader = SamSource(cram)
      val records = reader.toSeq
      records.size shouldBe 2
      records.head.name shouldBe "q1"
      reader.close()
    } finally {
      restoreGlobalState(snapshot)
    }
  }

  it should "allow explicit ref to override cramRefFasta when reading CRAM files" in {
    val snapshot = captureGlobalState()
    try {
      // Create a small reference sequence (500bp chr1)
      val refSequence = "N" * 500
      val refLen = refSequence.length

      // Create synthetic reference FASTA file
      val ref = makeTempFile("reference.", ".fasta")
      val refWriter = Io.toWriter(ref)
      refWriter.write(s">chr1\n")
      refWriter.write(refSequence + "\n")
      refWriter.close()

      // Create .fai index file for the reference
      val fai = Paths.get(ref.toString + ".fai")
      val faiWriter = Io.toWriter(fai)
      faiWriter.write(s"chr1\t$refLen\t6\t$refLen\t${refLen + 1}\n")
      faiWriter.close()

      // Set common args with a DIFFERENT path (which would fail if used)
      FgBioCommonArgs.args = new FgBioCommonArgs(cramRefFasta = Some(Paths.get("/nonexistent/ref.fasta")))

      // Create test data with custom sequence dictionary
      val dict = SequenceDictionary(SequenceMetadata(name="chr1", length=refLen))
      val builder = new SamBuilder(sd = Some(dict))
      builder.addPair(name="q1", start1=100, start2=200)
      val cram = makeTempFile("test.", ".cram")

      // Write to CRAM with explicit ref (overriding the broken common arg)
      val writer = SamWriter(cram, builder.header, ref = Some(ref))
      writer ++= builder.iterator
      writer.close()

      // Read from CRAM with explicit ref (overriding the broken common arg)
      val reader = SamSource(cram, ref = Some(ref))
      val records = reader.toSeq
      records.size shouldBe 2
      reader.close()
    } finally {
      restoreGlobalState(snapshot)
    }
  }
}
