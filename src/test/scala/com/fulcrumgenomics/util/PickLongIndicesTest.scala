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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.sopt.cmdline.ValidationException

class PickLongIndicesTest extends UnitSpec {

  "PickLongIndices.homopolymerLengthOk" should "correctly identify whether or not a sequence is ok based on homopolymer length" in {
    val picker = new PickLongIndices(length=12, numberOfIndices=100, output=Io.StdOut)

    // Check at the start of the sequence
    picker.homopolymerLengthOk("AAACCCGGGTTT".getBytes, longest=3) shouldBe true
    picker.homopolymerLengthOk("AAAACCCGGGTT".getBytes, longest=3) shouldBe false
    picker.homopolymerLengthOk("AAAACCCGGGTT".getBytes, longest=4) shouldBe true
    picker.homopolymerLengthOk("AAAAACCGGGTT".getBytes, longest=4) shouldBe false

    // Check at the end of the sequence
    picker.homopolymerLengthOk("AAACCCGGGTTT".getBytes, longest=3) shouldBe true
    picker.homopolymerLengthOk("AAACCCGGTTTT".getBytes, longest=3) shouldBe false
    picker.homopolymerLengthOk("AAACCCGGTTTT".getBytes, longest=4) shouldBe true
    picker.homopolymerLengthOk("AAACCCGTTTTT".getBytes, longest=4) shouldBe false

    // Check internal to the sequence
    picker.homopolymerLengthOk("AAACCCGGGTTT".getBytes, longest=3) shouldBe true
    picker.homopolymerLengthOk("AAACCCCGGTTT".getBytes, longest=3) shouldBe false
    picker.homopolymerLengthOk("AAACCCCGGTTT".getBytes, longest=4) shouldBe true
    picker.homopolymerLengthOk("AAACCCCCGTTT".getBytes, longest=4) shouldBe false

    // Check a bunch of random conditions that should pass
    picker.homopolymerLengthOk("AAACAAACAAAC".getBytes, longest=3) shouldBe true
    picker.homopolymerLengthOk("AAACAAACCAAA".getBytes, longest=3) shouldBe true
  }

  "PickLongIndices.isPalindrome" should "correctly identify palindrome vs. non-palindrome" in {
    val picker = new PickLongIndices(length = 12, numberOfIndices = 100, output = Io.StdOut)
    // Even length sequences
    picker.isPalindrome("AAAAAAAAAA".getBytes) shouldBe false
    picker.isPalindrome("AAAAATTTTT".getBytes) shouldBe true
    picker.isPalindrome("GGCCTTCCGG".getBytes) shouldBe false
    picker.isPalindrome("GGCCTAGGCC".getBytes) shouldBe true

    // Odd length sequences
    picker.isPalindrome("AAAAAAAAAAA".getBytes) shouldBe false
    picker.isPalindrome("AAAAAGTTTTT".getBytes) shouldBe true // middle base doesn't count
  }

  "PickLongIndices.AdapterRegex" should "correctly allow indices with a single block of Ns and not others" in {
    val picker = new PickLongIndices(length = 12, numberOfIndices = 100, output = Io.StdOut)
    picker.adapterRegex.matcher("ACGTNNNNACGT").matches shouldBe true
    picker.adapterRegex.matcher("ACGTNNNNNNNN").matches shouldBe true
    picker.adapterRegex.matcher("NNNNNNNNACGT").matches shouldBe true
    picker.adapterRegex.matcher("NNNNNNNNNNNN").matches shouldBe true

    picker.adapterRegex.matcher("ACGTTGATACGC").matches shouldBe false
    picker.adapterRegex.matcher("NNNNTTGANNNN").matches shouldBe false
  }

  "PickLongIndices" should "pick some indices given easy inputs" in {
    val picker = new PickLongIndices(length = 8, numberOfIndices = 96, editDistance = 3, output = Io.StdOut)
    val indices = picker.pickIndices
    indices should have size 96
    indices.foreach { index => index.length shouldBe 8 }
  }

  it should "pick some indices given easy inputs and some existing indices" in {
    val existing = Seq("ACGTACGT", "CAGTCAGT", "ATATATAT", "TATATATA", "ACACGCGC")
    val existingPath = makeTempFile("existing.", ".txt")
    Io.writeLines(existingPath, existing)
    val picker = new PickLongIndices(length = 8, numberOfIndices = 96, editDistance = 3, output = Io.StdOut, existing=Some(existingPath))
    val indices = picker.pickIndices
    indices should have size 96
    indices.foreach { index => index.length shouldBe 8 }

    val set = indices.map(new String(_)).toSet
    existing.foreach(e => set.contains(e) shouldBe true)
  }

  it should "run end to end and produce an output file" in {
    val out = makeTempFile("indices.", ".txt")
    val picker = new PickLongIndices(length = 8, numberOfIndices = 96, editDistance = 3, output = out)
    picker.execute()
    out.toFile.exists() shouldBe true
    val lines = Io.readLines(out).toSeq
    lines should have size 97 // 96 plus header
  }

  it should "reject invalid option combinations" in {
    an[ValidationException] should be thrownBy { new PickLongIndices(length=0, numberOfIndices=24, editDistance=3,   output=Io.StdOut) }
    an[ValidationException] should be thrownBy { new PickLongIndices(length=8, numberOfIndices=0 , editDistance=3,   output=Io.StdOut) }
    an[ValidationException] should be thrownBy { new PickLongIndices(length=8, numberOfIndices=24, editDistance= -1, output=Io.StdOut) }
    an[ValidationException] should be thrownBy { new PickLongIndices(length=8, numberOfIndices=24, editDistance=9, maxHomopolymer = -1, output=Io.StdOut) }

    val existingPath = makeTempFile("existing.", ".txt")
    val existing = Seq("AAAAA", "AAAAATTTTT")
    Io.writeLines(existingPath, existing)
    an[ValidationException] should be thrownBy { new PickLongIndices(length=8, numberOfIndices=24, editDistance=3, existing=Some(existingPath), output=Io.StdOut) }
  }
}
