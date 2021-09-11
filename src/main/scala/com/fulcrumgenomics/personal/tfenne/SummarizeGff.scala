/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.personal.tfenne

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.{LazyLogging, SimpleCounter}
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.GeneAnnotations.GeneBiotype
import com.fulcrumgenomics.util.NcbiRefSeqGffSource
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

@clp(group=ClpGroups.Personal, description=
  """
    |Takes in a RefSeq GFF and produces some summary statistics.
  """)
class SummarizeGff(@arg(flag='i', doc="Input RefSeq GFF") input: FilePath,
                   @arg(flag='r', doc="Reference sequence fasta.") ref: PathToFasta,
                   @arg(flag='o', doc="Output summary file.") output: FilePath) extends FgBioTool with LazyLogging {

  override def execute(): Unit = {
    val dict = SequenceDictionary.extract(ref)
    val parser = NcbiRefSeqGffSource(input, includeXs=true, dict)

    val singleExonCounts = new SimpleCounter[String]
    val multiExonCounts  = new SimpleCounter[String]

    parser.foreach { gene =>
      val biotype = gene.biotype.map(_.toString).getOrElse("unknown")
      val multi = gene.loci.exists(_.transcripts.exists(_.exons.size > 1))
      if (multi) multiExonCounts.count(biotype) else singleExonCounts.count(biotype)
    }

    println("biotype\tsingle_exon_genes\tmulti_exon_genes")
    GeneBiotype.values.map(_.toString).sorted.foreach { biotype =>
      println(s"${biotype}\t${singleExonCounts.countOf(biotype)}\t${multiExonCounts.countOf(biotype)}")
    }
  }
}
