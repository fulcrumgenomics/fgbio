/*
 * The MIT License
 *
 * Copyright (c) 2022 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.fasta

import com.fulcrumgenomics.FgBioDef.PathToSequenceDictionary
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.CommonsDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt._
import com.fulcrumgenomics.util.Io

import scala.collection.immutable.IndexedSeq
import scala.collection.mutable.{ListBuffer, Builder}

@clp(description =
  """
    |Sorts a sequence dictionary file in the order of another sequence dictionary.
    |
    |The inputs are to two `*.dict` files.  One to be sorted, and the other to provide the order for the sorting.
    |
    |If there is a contig in the input dictionary that is not in the sorting dictionary, that contig will be appended
    |to the end of the sequence dictionary, unless appendMissingContigs is set to false, in which case such contigs will be ignored.
    |
    |If there is a contig in the sorting dictionary that is not in the input dictionary, that contig will be ignored.
    |
    |The output will be a "sequence dictionary", which is a valid SAM file, containing the version header line and one
    |line per contig.  The fields of the entries in this dictionary will be the same as in input, but in the order of
    |`sortDictionary`.
  """,
  group = ClpGroups.Fasta)
class SortSequenceDictionary
(@arg(flag='i', doc="Input sequence dictionary file to be sorted.") val input: PathToSequenceDictionary,
 @arg(flag='d', doc="Input sequence dictionary file containing contigs in the desired sort order.") val sortDictionary: PathToSequenceDictionary,
 @arg(flag='o', doc="Output sequence dictionary file.") val output: PathToSequenceDictionary,
 @arg(doc="Append input contigs that have no matching contig in the sort dictionary to the end of the output dictionary. If false, such contigs will not be output.") val appendMissingContigs: Boolean = true,
) extends FgBioTool with LazyLogging {
  Io.assertReadable(input)
  Io.assertReadable(sortDictionary)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
      val inputDict = SequenceDictionary(input)
      val sortOrderDict = SequenceDictionary(sortDictionary)

      val newOrder = IndexedSeq.newBuilder[SequenceMetadata]
      sortOrderDict.foreach { sortMeta =>
          sortMeta.allNames.find { name => inputDict.contains(name) } match {
            case Some(name) => {
              val inMeta = inputDict.get(name).getOrElse(throw new java.lang.RuntimeException(s"Contig '${sortMeta.name}' alias '${name}' not found in input dictionary after initially being found (This should never happen)"))
              newOrder += inMeta.copy()
            }
            case None => logger.info(s"Contig '${sortMeta.name}' corresponded to no contig in input dictionary, skipping")
          }
      }

      val tmpDict = {
          val metadata = newOrder.result().zipWithIndex.map {
            case (meta, index) => meta.copy(index=index)
          }.toSeq
          SequenceDictionary(metadata:_*)
      }

      inputDict.foreach { inMeta =>
        if (!tmpDict.contains(inMeta.name)) {
          val skipBehavior = if (appendMissingContigs) "appending." else "skipping."
          logger.info(s"Contig '${inMeta.name}' was not found in sort order dictionary: ${skipBehavior}")
          if (appendMissingContigs) {
              newOrder += inMeta.copy()
          }
        }
      }
      val metadata = newOrder.result().zipWithIndex.map {
          case (meta, index) => meta.copy(index=index)
      }.toSeq
      SequenceDictionary(metadata:_*).write(output)
  }
}
