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
 */

package com.fulcrumgenomics.bam.api

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.SamSource.{DefaultUseAsyncIo, DefaultValidationStringency}
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Io
import htsjdk.samtools.GenomicIndexUtil

import java.nio.file.Files
import java.util.concurrent.{Callable, Executors, TimeUnit}
import scala.util.Random

class SamIoTest extends UnitSpec {
  "SamWriter/SamSource" should "write some records to a BAM on disk and read them back" in {
    val builder = new SamBuilder()
    builder.addPair(name="q1", start1=100, start2=300)
    builder.addPair(name="q4", start1=200, start2=400)
    builder.addPair(name="q3", start1=300, start2=500)
    builder.addPair(name="q2", start1=400, start2=600)

    val bam = makeTempFile("test.", ".bam")
    val out = SamWriter(bam, builder.header, sort=Some(SamOrder.Queryname))
    out ++= builder.iterator
    out.close()

    val in = SamSource(bam)
    in.toSeq.map(_.name) should contain theSameElementsInOrderAs  Seq("q1","q1","q2","q2","q3","q3","q4","q4")
  }

  it should "write some records to a SAM text file on disk and read them back" in {
    val builder = new SamBuilder()
    builder.addPair(name="q1", start1=100, start2=300)
    builder.addPair(name="q4", start1=200, start2=400)
    builder.addPair(name="q3", start1=300, start2=500)
    builder.addPair(name="q2", start1=400, start2=600)

    val sam = makeTempFile("test.", ".sam")
    val out = SamWriter(sam, builder.header, sort=Some(SamOrder.Coordinate))
    out ++= builder.iterator
    out.close()

    val in = SamSource(sam)
    in.toSeq.map(_.start) should contain theSameElementsInOrderAs Seq(100, 200, 300, 300, 400, 400, 500, 600)

    Io.readLines(sam).next() should startWith ("@HD")
  }

  it should "sort to disk and read records back" in {
    val builder = new SamBuilder(sort=None)
    val random  = new Random(42)
    Range.inclusive(1, 1000).foreach { _ =>
      builder.addPair(contig=random.nextInt(23), start1=random.nextInt(1000000), start2=random.nextInt(1000000))
    }

    val tmp = makeTempFile("sorted.", ".bam")
    val out = SamWriter(tmp, builder.header, sort=Some(SamOrder.Coordinate), maxRecordsInRam=500)
    out ++= builder.iterator
    out.close()
    val recs = SamSource(tmp).iterator.toSeq
    recs.size shouldBe builder.size
    recs.grouped(2).foreach { case Seq(r1, r2) => r1.refIndex < r2.refIndex || r1.start <= r2.start shouldBe true }
  }

  it should "support streaming through named pipes" in {
    // Make a named pipe
    val pipe = makeTempFile("pipey.", ".mcpiperton")
    Files.delete(pipe)
    val mkfifo = new ProcessBuilder("mkfifo", pipe.toAbsolutePath.toString).start()
    mkfifo.waitFor(2, TimeUnit.SECONDS) shouldBe true
    mkfifo.exitValue() shouldBe 0

    val builder = new SamBuilder(readLength=100)
    val _ = forloop (from=1, until=5001) { i => val _ = builder.addFrag(name=s"q$i", contig=0, start=i) }

    // Figure up an executor to pump data into the pipe
    val exec = Executors.newSingleThreadExecutor()
    val future = exec.submit(new Callable[Int] {
      override def call(): Int = {
        val out = SamWriter(pipe, header=builder.header)
        out ++= builder
        out.close()
        builder.size
      }
    })

    // Then try reading from the pipe
    val in  = SamSource(pipe)
    var count = 0
    in.foreach { rec =>
      rec.name shouldBe s"q${rec.start}"
      count += 1
    }

    in.safelyClose()
    count shouldBe future.get()
  }

  it should "not explode when asked to index a file with chroms too long" in {
    val dict    = SequenceDictionary(SequenceMetadata(name="chr1", length=GenomicIndexUtil.BIN_GENOMIC_SPAN + 100000))
    val builder = new SamBuilder(sd=Some(dict), sort=Some(SamOrder.Coordinate))
    builder.addFrag(start=GenomicIndexUtil.BIN_GENOMIC_SPAN + 500)

    val bam = makeTempFile("long.", ".bam")
    val out = SamWriter(bam, header=builder.header, index=true)
    out ++= builder
    out.close()

    val in = SamSource(bam)
    in.indexed shouldBe false
  }

  "HeaderHelper" should "provide tidy access to things in the header" in {
    val builder = new SamBuilder()
    val source = builder.toSource
    source.readGroups should contain theSameElementsInOrderAs source.header.getReadGroups.toSeq
    source.dict shouldBe source.dict
    source.programGroups should contain theSameElementsInOrderAs source.header.getProgramRecords.toSeq
  }

  "SamSource" should "iterate lazily and not suck all records into memory when using filter() or map()" in {
    val builder = new SamBuilder(readLength=10, baseQuality=20)
    Range(0, 10).foreach { _ => builder.addFrag(start=100) }
    val source = builder.toSource

    var filterCount: Int = 0
    var mapCount: Int = 0

    val xs = source.filter { _ =>
      filterCount += 1
      true
    }.map { r =>
      mapCount += 1
      r
    }

    filterCount shouldBe 0
    mapCount shouldBe 0
    xs.toList
    filterCount shouldBe 10
    mapCount shouldBe 10
  }

  it should "allow reading from a stream of SAM bytes" in {
    val builder = new SamBuilder()
    builder.addPair(name = "q1", start1 = 100, start2 = 300)
    builder.addPair(name = "q4", start1 = 200, start2 = 400)
    builder.addPair(name = "q3", start1 = 300, start2 = 500)
    builder.addPair(name = "q2", start1 = 400, start2 = 600)

    val sam = makeTempFile(getClass.getSimpleName, ".sam")
    val out = SamWriter(sam, builder.header, sort = Some(SamOrder.Coordinate))
    builder.foreach(out.write)
    out.close()

    val source = SamSource(
      stream     = Io.toInputStream(sam),
      ref        = None,
      async      = DefaultUseAsyncIo,
      stringency = DefaultValidationStringency,
      factory    = SamRecord.Factory
    )

    source.indexed shouldBe false
    source.toSeq.map(_.start) should contain theSameElementsInOrderAs Seq(100, 200, 300, 300, 400, 400, 500, 600)
    source.close()
  }

  "SamSource.query" should "query a BAM, accounting for queries before and after the contig start and end" in {
      val queryType = QueryType.Overlapping
      val builder = new SamBuilder(readLength=10, baseQuality=20, sort=Some(SamOrder.Coordinate))
      Range(0, 10).foreach { _ => builder.addFrag(start=100) }
      val source = builder.toSource

      // test a query before the contig start
      source.query("chr1", 0, 1000, queryType).length shouldBe 10
      source.query("chr1", -100, 1000, queryType).length shouldBe 10

      // test a query after the contig end
      val contigEnd = builder.dict("chr1").length
      source.query("chr1", 100, contigEnd+1, queryType).length shouldBe 10
      source.query("chr1", 100, contigEnd+100, queryType).length shouldBe 10

      // test a query both before and after the contig start and end respectively
      source.query("chr1", -100, contigEnd+100, queryType).length shouldBe 10

      // at the start and end     
      source.query("chr1", 1, 1000, queryType).length shouldBe 10 
      source.query("chr1", 100, contigEnd, queryType).length shouldBe 10
      source.query("chr1", 1, contigEnd, queryType).length shouldBe 10

      // exception when the contig does not exist
      an[NoSuchElementException] should be thrownBy source.query("contig-does-not-exist", 1, 1000, queryType)
  }
}
