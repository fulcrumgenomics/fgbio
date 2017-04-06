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
 *
 */

package com.fulcrumgenomics.basecalling

import com.fulcrumgenomics.FgBioDef.FilePath
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.illumina.RunInfo
import com.fulcrumgenomics.util._
import dagr.sopt.{arg, clp}

@clp(group=ClpGroups.Basecalling, description=
  """
    |Extracts information about an Illumina sequencing run from the RunInfo.xml.
    |
    |The output file will contain a header column and a single column containing the following rows:
    |1. run_barcode: the unique identifier for the sequencing run and flowcell, stored as "<instrument-name>_<flowcell-barcode>".
    |2. flowcell_barcode: the flowcell barcode.
    |3. instrument_name: the instrument name.
    |4. run_date: the date of the sequencing run.
    |5. read_structure: the description of the logical structure of cycles within the sequencing run, including which cycles
    |   correspond to sample barcodes, molecular barcodes, template bases, and bases that should be skipped.
    |6. number_of_lanes: the number of lanes in the flowcell.
  """)
class ExtractIlluminaRunInfo
(
  @arg(flag="i", doc="The input RunInfo.xml typically found in the run folder.") val input: FilePath,
  @arg(flag="o", doc="The output file.") val output: FilePath
) extends FgBioTool {

  Io.assertReadable(input)
  Io.assertCanWriteFile(output)

  override def execute(): Unit = {
    Metric.write(output, RunInfo(runInfo=input))
  }
}
