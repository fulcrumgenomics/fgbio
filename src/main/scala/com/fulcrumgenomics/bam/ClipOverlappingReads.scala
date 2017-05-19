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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.ClpGroups
import com.fulcrumgenomics.sopt.{arg, clp}

@deprecated(since="0.2.0", message="Use ClipBam instead")
@clp(group = ClpGroups.SamOrBam, description=
  """
    |Clips reads from the same template to eliminate overlap between the reads. Ensures that downstream
    |processes, particularly variant calling, cannot double-count evidence from the same template when
    |both reads span a variant site in the same template.
    |
    |Clipping is only performed on `FR` read pairs, and is implemented by clipping approximately half the
    |overlapping bases from each read.  By default hard clipping is performed; soft-clipping may be
    |substituted using the `--soft-clip` parameter.
    |
    |Secondary alignments and supplemental alignments are not clipped, but are passed through into the
    |output.
    |
    |If the input BAM is neither `queryname` sorted nor `query` grouped, it will be sorted into queryname
    |order so that clipping can be performed on both ends of a pair simultaneously and so that mate
    |pair information can be reset across all reads for the template.  Post-clipping the reads are
    |resorted into coordinate order, any existing `NM`, `UQ` and `MD` tags are repaired, and the output is
    |written in coordinate order.
  """)
class ClipOverlappingReads
( @arg(flag='i', doc="Input SAM or BAM file of aligned reads in coordinate order.") input: PathToBam,
  @arg(flag='o', doc="Output SAM or BAM file.") output: PathToBam,
  @arg(flag='r', doc="Reference sequence fasta file.") ref: PathToFasta,
  @arg(flag='s', doc="Soft clip reads instead of hard clipping.") softClip: Boolean = false,
  @arg(flag='a', doc="Automatically clip extended attributes that are the same length as bases.") autoClipAttributes: Boolean = false
) extends ClipBam(input=input, output=output, ref=ref, softClip=softClip, autoClipAttributes=autoClipAttributes)
