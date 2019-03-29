/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.fastq

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.{Io, Sorter}

/** Companion object for the SortFastq tool. */
object SortFastq {
  private class FastqCodec extends Sorter.Codec[FastqRecord] {
    /** Encodes the record directly to fastq. */
    override def encode(rec: FastqRecord): Array[Byte] = rec.toString.getBytes

    /** Decodes the record from the fastq string. */
    override def decode(bs: Array[Byte], start: Int, length: Int): FastqRecord = {
      val lines = new String(bs, start, length).lines
      val name       = lines.next().substring(1)
      val bases      = lines.next()
      val qualHeader = lines.next()
      val quals      = lines.next()

      FastqRecord(name=name, bases=bases, quals=quals)
    }
  }

  private class SortKey(val name: String, val number: Byte) extends Ordered[SortKey] {
    override def compare(that: SortKey): Int = {
      var retval = this.name.compareTo(that.name)
      if (retval == 0) retval = this.number - that.number
      retval
    }
  }
}

@clp(group=ClpGroups.Personal, description=
  """
    |Sorts a FASTQ file.  Sorts the records in a FASTQ file based on the lexicographic ordering
    |of their read names.  Input and output files can be either uncompressed or gzip-compressed.
  """)
class SortFastq
( @arg(flag='i', doc="Input fastq file.") input: PathToFastq,
  @arg(flag='o', doc="Output fastq file.") output: PathToFastq,
  @arg(flag='m', doc="Maximum records to keep in RAM at one time.") val maxRecordsInRam: Int = 5e5.toInt
) extends FgBioTool {
  import SortFastq._

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    val in = FastqSource(input)
    val sorter = new Sorter[FastqRecord, SortKey](
      maxObjectsInRam = maxRecordsInRam,
      codec           = new FastqCodec,
      keyfunc         = r => new SortKey(r.name, r.readNumber.getOrElse(1).toByte)
    )

    sorter ++= in
    in.safelyClose()

    val out = FastqWriter(output)
    out ++= sorter.iterator
    out.close()
  }
}
