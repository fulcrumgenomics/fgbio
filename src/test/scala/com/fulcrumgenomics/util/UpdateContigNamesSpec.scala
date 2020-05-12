/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef.PathToSequenceDictionary
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.testing.UnitSpec

import scala.collection.mutable.ListBuffer

trait UpdateContigNamesSpec extends UnitSpec {

  protected def toSequenceMetadata(name: String, alias: String*): SequenceMetadata = {
    SequenceMetadata(name=name, length=100000000, aliases=alias)
  }

  protected def dict(skipLast: Boolean = false): SequenceDictionary = {
    val infos = ListBuffer[SequenceMetadata](
      toSequenceMetadata(name="chr1", "NC_000001.10"),
      toSequenceMetadata(name="chr2", "NC_000002.10"),
      toSequenceMetadata(name="chr3", "NC_000003.10")
    )
    if (!skipLast) infos += toSequenceMetadata(name="chr4", "NC_000004.10")
    SequenceDictionary(infos.toSeq:_*)
  }

  protected def pathToSequenceDictionary(skipLast: Boolean = false): PathToSequenceDictionary = {
    val path = makeTempFile("test.", "in.dict")
    dict(skipLast=skipLast).write(path)
    path
  }
}
