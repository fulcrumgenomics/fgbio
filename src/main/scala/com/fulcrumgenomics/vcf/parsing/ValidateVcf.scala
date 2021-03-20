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

package com.fulcrumgenomics.vcf.parsing

import com.fulcrumgenomics.FgBioDef.{BetterBufferedIteratorScalaWrapper, PathToVcf}
import com.fulcrumgenomics.cmdline.{ClpGroups, FgBioTool}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.sopt.{arg, clp}
import com.fulcrumgenomics.util.Io
import com.fulcrumgenomics.vcf.api.{Variant, VcfHeader, VcfHeaderEntry}
import com.fulcrumgenomics.vcf.parsing.header.VcfHeaderBuilder
import com.fulcrumgenomics.vcf.parsing.header.VcfHeaderParser._
import com.fulcrumgenomics.vcf.parsing.util.ParseResult.LogFailure
import com.fulcrumgenomics.vcf.parsing.variant.VariantParser
import fastparse.Parsed.{Failure, Success}
import fastparse.internal.{Logger => FastParseLogger}
import fastparse.{Parsed, parse}

@clp(
  description =
    """
      |Validates a VCF by trying to parse it.
      |
      |Currently only supports VCF specification 4.3.
      |
      |Many things are NOT validated, user beware.
      |
    """,
  group=ClpGroups.VcfOrBcf
)
class ValidateVcf
( @arg(flag='i', doc="The input VCF.") val input: PathToVcf,
) extends FgBioTool with LazyLogging {

  private val logBuffer = collection.mutable.Buffer.empty[String]
  implicit private val fastParseLogger: FastParseLogger = FastParseLogger(logBuffer.append)

  override def execute(): Unit = {
    val headerBuilder = new VcfHeaderBuilder()
    val lineIter      = Io.readLines(path=input).bufferBetter
    var lineNumber    = 1

    // Parse the header lines
    lineIter.takeWhile(_.startsWith("##")).foreach { line=>
      handleResult[VcfHeaderEntry](
        result     = parse(line, headerBuilder.parse(_)),
        lineNumber = lineNumber
      )
      lineNumber += 1
    }

    // Parse the samples
    val sampleNames = lineIter.nextOption() match {
      case None       => fail("Pre-mature end of file: expected a #CHROM line")
      case Some(line) => handleResult(result=parse(line, samples(_)), lineNumber=lineNumber)
    }
    lineNumber += 1

    // Build the header
    val header: VcfHeader = headerBuilder.header(samples=sampleNames)

    // Parse the variants
    val parser: VariantParser = new VariantParser(header=header)
    lineIter.foreach { line =>
      handleResult[Variant](
        result     = parse(line, parser.variant(_)),
        lineNumber = lineNumber
      )
      lineNumber += 1
    }
  }

  private def handleResult[T](result: Parsed[T], lineNumber: Int): T = {
    result match {
      case Success(entry, _) => entry
      case failure: Failure  =>
        // Log an error
        val message = logger.fail(failure=failure, lineNumber=lineNumber, logBuffer=logBuffer.iterator)
        logBuffer.clear()
        fail(message)
    }
  }
}
