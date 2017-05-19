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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.util.Sorter.Codec

private class StringCodec extends Codec[String] {
  override def encode(a: String): Array[Byte] = a.getBytes
  override def decode(bs: Array[Byte], start: Int, length: Int): String = new String(bs, start, length)
}

private case class OrderedString(s: String) extends Ordered[OrderedString] {
  override def compare(that: OrderedString): Int = s.compare(that.s)
}


class SorterTest extends UnitSpec {
  "Sorter" should "sort a bunch of Strings" in {
    val sorter = new Sorter(10, new StringCodec, OrderedString.apply)
    val strings = Seq("hello", "Hello", "hey!", "whoa", "hola", "goodbye", "seeya", "G'day", "bienvenue", "guttentarg", "adios")
    sorter ++= strings
    val sorted = sorter.iterator.toSeq
    sorted shouldBe strings.sorted
  }

  it should "work if no items are added" in {
    val sorter = new Sorter(10, new StringCodec, OrderedString.apply)
    sorter.iterator.toSeq shouldBe Seq.empty
  }
}
