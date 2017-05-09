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
package com.fulcrumgenomics.cmdline

import com.fulcrumgenomics.sopt.cmdline.ClpGroup

/** Groups for organizing command line programs for display. */
object ClpGroups {
  class _PersonalGroup extends ClpGroup {
    override val name: String = "Personal"
    override val description: String =  "Various personal programs (not supported)."
    override val rank = Integer.MAX_VALUE
  }

  class _SamOrBamGroup extends ClpGroup {
    override val name: String = "SAM/BAM"
    override val description: String = "Tools for manipulating SAM, BAM, or related data."
  }

  class _FastaGroup extends ClpGroup {
    override val name: String = "FASTA"
    override val description: String = "Tools for manipulating FASTA files."
  }

  class _FastqGroup extends ClpGroup {
    override val name: String = "FASTQ"
    override val description: String = "Tools for manipulating FASTQ files."
  }

  class _UtilitiesGroup extends ClpGroup {
    override val name: String = "Utilities"
    override val description: String = "Various utility programs."
  }

  class _Umi extends ClpGroup {
    override def name: String = "Unique Molecular Identifiers (UMIs)"
    override def description: String = "Tools for manipulating UMIs & reads tagged with UMIs"
  }

  class _VcfOrBcfGroup extends ClpGroup {
    override val name: String = "VCF/BCF"
    override val description: String = "Tools for manipulating VCF, BCF, or related data."
  }

  class _Basecalling extends ClpGroup {
    override val name: String = "Basecalling"
    override val description: String = "Tools for manipulating basecalling data."
  }

  final val Personal    = classOf[_PersonalGroup]
  final val SamOrBam    = classOf[_SamOrBamGroup]
  final val Fasta       = classOf[_FastaGroup]
  final val Fastq       = classOf[_FastqGroup]
  final val Umi         = classOf[_Umi]
  final val Utilities   = classOf[_UtilitiesGroup]
  final val VcfOrBcf    = classOf[_VcfOrBcfGroup]
  final val Basecalling = classOf[_Basecalling]
}

