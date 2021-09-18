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

package com.fulcrumgenomics.vcf.parsing.header

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import enumeratum.EnumEntry

import scala.collection.immutable

/** Trait for common keys in a meta header VCF entry/line. */
sealed trait MetaHeaderKey extends EnumEntry

object MetaHeaderKey extends FgBioEnum[MetaHeaderKey] {
  case object ID extends MetaHeaderKey
  case object Number extends MetaHeaderKey
  case object Type extends MetaHeaderKey
  case object Description extends MetaHeaderKey
  case object Source extends MetaHeaderKey
  case object Version extends MetaHeaderKey
  override def values: immutable.IndexedSeq[MetaHeaderKey] = findValues

  implicit def metaHeaderKeyToString(key: MetaHeaderKey): String = key.entryName
}