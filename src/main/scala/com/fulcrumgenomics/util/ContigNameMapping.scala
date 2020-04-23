/*
 * The MIT License
 *
 * Copyright (c) 2020 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.FgBioDef._

object  ContigNameMapping {
    /**
      * Parses a contig name mappings file that is a tab seperated source -> target mapping. EX:
      *
      * """
      * 1   chr1
      * 2   chr2    Chr2
      * 3   chr3
      * 3   Chr3
      * """
      *
      * Multiple targets can be specified on one line, or on multiple lines.
      */
    def parse(nameMapping: FilePath): Map[String, Seq[String]] = {
        Io.readLines(nameMapping).map { line =>
            line.split('\t').toList match {
                case head::tail if tail.length > 0 => (head, tail)
                case _ => throw new IllegalArgumentException(s"Malformed line: expected at least two columns: $line")
            }
        }.toSeq.groupBy(_._1).map(k => k._1 -> k._2.map(_._2).flatten.distinct.toList)
    }
}