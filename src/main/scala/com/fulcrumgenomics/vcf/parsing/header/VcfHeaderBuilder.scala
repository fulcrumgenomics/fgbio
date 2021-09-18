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

import com.fulcrumgenomics.vcf.api._
import com.fulcrumgenomics.vcf.parsing.util.ParseResult.success
import fastparse.NoWhitespace._
import fastparse.{P, Pass}

import scala.collection.mutable.ArrayBuffer


/** Buffers VCF lines of the same type.
  *
  * @param lineName the name of the line type (eg. contig, INFO, FORMAT, FILTER)
  * @param keySource the source of key used to uniquely identify a line (eg. `name` for contig, `ID` for INFO, ...)
  * @tparam LineType the type of [[VcfHeaderEntry]]
  */
private class VcfLineBuffer[LineType<:VcfKeyValHeader]
(lineName: String, keySource: String) extends Iterable[LineType] {
  import com.fulcrumgenomics.vcf.parsing.util.ParseResult.fail

  private val buffer = ArrayBuffer[LineType]()

  /** Parses a given line.  Fails if the a previous line has the same identifier. */
  def parse[_: P](line: LineType): P[LineType] = {
    val id: String = line.id
    buffer.find(_.id == id) match {
      case None    => this.buffer += line; success(line)
      case Some(_) => fail(startIndex=0, f"duplicate $lineName with $keySource `$id`")
    }
  }

  /** Returns an interator over the buffered lines. */
  def iterator: Iterator[LineType] = this.buffer.iterator
}

/** Helps build a [[VcfHeader]] entry/line-by-entry. */
class VcfHeaderBuilder() {
  private val contigs = new VcfLineBuffer[VcfContigHeader]         (lineName="contig", keySource="name")
  private val infos   = new VcfLineBuffer[VcfInfoHeader]           (lineName="INFO",   keySource="id")
  private val formats = new VcfLineBuffer[VcfFormatHeader]         (lineName="FORMAT", keySource="id")
  private val filters = new VcfLineBuffer[VcfFilterHeader]         (lineName="FILTER", keySource="id")
  private val alts    = new VcfLineBuffer[VcfAlternateAlleleHeader](lineName="ALT", keySource="id")
  private val others  = new VcfLineBuffer[VcfGeneralHeader]        (lineName="other",  keySource="id")

  import com.fulcrumgenomics.vcf.parsing.header.VcfHeaderParser.{header => parseHeader}

  /** Parse a single header entry/line. */
  def parse[_: P]: P[VcfHeaderEntry] = parseHeader.flatMap {
    case contig: VcfContigHeader       => contigs.parse(contig)
    case info: VcfInfoHeader           => infos.parse(info)
    case format: VcfFormatHeader       => formats.parse(format)
    case filter: VcfFilterHeader       => filters.parse(filter)
    case alt: VcfAlternateAlleleHeader => alts.parse(alt)
    case other: VcfGeneralHeader       => others.parse(other)
  }

  /** Build a VCF header with the given sample names. */
  def header(samples: Seq[String]): VcfHeader = {
    VcfHeader(
      contigs = contigs.zipWithIndex.map { case (contig, index) => contig.copy(index=index) }.toIndexedSeq,
      infos   = infos.toIndexedSeq,
      formats = formats.toIndexedSeq,
      filters = filters.toIndexedSeq,
      alts    = alts.toIndexedSeq,
      others  = others.toIndexedSeq,
      samples = samples.toIndexedSeq
    )
  }
}