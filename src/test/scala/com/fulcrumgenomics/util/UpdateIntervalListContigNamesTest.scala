/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

import com.fulcrumgenomics.FgBioDef.PathToIntervals
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import htsjdk.samtools.util.{Interval, IntervalList}

import scala.collection.mutable.ListBuffer

class UpdateIntervalListContigNamesTest extends UpdateContigNamesSpec {

  private def inIntervalList: IntervalList = {
    import com.fulcrumgenomics.fasta.Converters.ToSAMSequenceDictionary
    val intervals = ListBuffer[Interval]()
    intervals += new Interval("NC_000001.10", 123, 500)
    intervals += new Interval("NC_000002.10", 444, 888, true, null)
    intervals += new Interval("NC_000003.10", 8284, 10000, true, "name")
    intervals += new Interval("NC_000004.10", 112345, 223456)

    val dict = SequenceDictionary(
      intervals.map(i => SequenceMetadata(name=i.getContig, length=100000000)).toSeq:_*
    )
    val intervalList = new IntervalList(dict.asSam)
    intervals.foreach(intervalList.add)
    intervalList
  }

  private def outIntervalList(skipLast: Boolean = false): Seq[Interval] = {
    val intervals = ListBuffer[Interval]()
    intervals += new Interval("chr1", 123, 500)
    intervals += new Interval("chr2", 444, 888, true, null)
    intervals += new Interval("chr3", 8284, 10000, true, "name")
    if (!skipLast) intervals += new Interval("chr4", 112345, 223456)
    intervals.toSeq
  }

  private val pathToInIntervals: PathToIntervals = {
    val path = makeTempFile("test.", "in.interval_list")
    inIntervalList.write(path)
    path
  }

  private def slurp(path: PathToIntervals): Seq[Interval] = {
    val source = IntervalListSource(path)
    val intervals = source.toIndexedSeq
    source.close()
    intervals
  }

  "UpdateIntervalListContigNames" should "update an Interval List" in {
    val output = makeTempFile("test.", "out.interval_list")
    val tool = new UpdateIntervalListContigNames(
      input  = pathToInIntervals,
      dict   = pathToSequenceDictionary(),
      output = output
    )

    executeFgbioTool(tool)

    slurp(output) should contain theSameElementsInOrderAs outIntervalList()
  }

  it should "throw an exception if there are missing source contigs" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateIntervalListContigNames(
      input  = pathToInIntervals,
      dict   = pathToSequenceDictionary(skipLast=true),
      output = output
    )

    val ex = intercept[Exception] { executeFgbioTool(tool) }
    ex.getMessage should include ("Did not find contig")
  }

  it should "skip missing source contigs when using --skip-missing" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateIntervalListContigNames(
      input       = pathToInIntervals,
      dict        = pathToSequenceDictionary(skipLast=true),
      output      = output,
      skipMissing = true
    )

    executeFgbioTool(tool)

    slurp(output) should contain theSameElementsInOrderAs outIntervalList(skipLast=true)
  }
}