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
    |to the end of the sequence dictionary in the same relative order to other appended contigs as in the input dictionary.
    |Missing contigs can be omitted by setting `--skip-missing-contigs` to true.
    |
    |If there is a contig in the sorting dictionary that is not in the input dictionary, that contig will be ignored.
    |
    |The output will be a sequence dictionary, containing the version header line and one
    |line per contig.  The fields of the entries in this dictionary will be the same as in input, but in the order of
    |`--sort-dictionary`.
  """,
  group = ClpGroups.Fasta)
class SortSequenceDictionary
(@arg(flag='i', doc="Input sequence dictionary file to be sorted.") val input: PathToSequenceDictionary,
 @arg(flag='d', doc="Input sequence dictionary file containing contigs in the desired sort order.") val sortDictionary: PathToSequenceDictionary,
 @arg(flag='o', doc="Output sequence dictionary file.") val output: PathToSequenceDictionary,
 @arg(doc="Skip input contigs that have no matching contig in the sort dictionary rather than appending to the end of the output dictionary.") val skipMissingContigs: Boolean = false,
) extends FgBioTool with LazyLogging {
  
  Io.assertReadable(input)
  Io.assertReadable(sortDictionary)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
      val inputDict     = SequenceDictionary(input)
      val sortOrderDict = SequenceDictionary(sortDictionary)

      // Iterate through the sort dictionary collecting metas from the input that match by name
      val metasBuilder = IndexedSeq.newBuilder[SequenceMetadata]
      sortOrderDict.foreach { sortMeta =>
          sortMeta.allNames.find { name => inputDict.contains(name) } match {
            case Some(name) =>  metasBuilder += inputDict(name)
            case None => logger.info(s"Contig '${sortMeta.name}' corresponded to no contig in input dictionary, skipping")
          }
      }

      // build a dictionary from the input contigs found in the sort dictionary
      val metasFoundInSortDictDict = {
          val metadata = metasBuilder.result().zipWithIndex.map {
            case (meta, index) => meta.copy(index=index)
          }.toSeq
          SequenceDictionary(metadata:_*)
      }

      // maybe append input contigs not found in the sort dictionary.  Their index will be reset after aggregation.
      inputDict.foreach { inMeta =>
        if (!metasFoundInSortDictDict.contains(inMeta.name)) {
          val skipBehavior = if (skipMissingContigs) "skipping." else "appending."
          logger.warning(s"Contig '${inMeta.name}' was not found in sort order dictionary: $skipBehavior")
          // Append if desired. The index will be reset later.
          if (!skipMissingContigs) {
              metasBuilder += inMeta.copy()
          }
        }
      }
      // Finally we have all the contigs, so reset the index and write out the dictionary. 
      val finalMetadataDict = metasBuilder.result().zipWithIndex.map {
          case (meta, index) => meta.copy(index=index)
      }.toSeq
      SequenceDictionary(finalMetadataDict:_*).write(output)
  }
}
