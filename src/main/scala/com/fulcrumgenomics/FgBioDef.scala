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

package com.fulcrumgenomics

import com.fulcrumgenomics.bam.api.SamSource
import com.fulcrumgenomics.commons.CommonsDef
import com.fulcrumgenomics.commons.io.PathUtil
import enumeratum.{Enum, EnumEntry}

/** FgBioDef that is just an extension of the commons CommonsDef. Here in case
  * there's a need for FgBio specific defs later.
  */
object FgBioDef extends CommonsDef {
  /** Extends this trait in your enumeration object to enable the case objects to be created on the command line.
    * You should implement the [[Enum#values]] method in your object, for example:
    * `def values: IndexedSeq[T] = findValues`.
    * */
  trait FgBioEnum[T<:EnumEntry] extends Enum[T] {
    // NB: we cannot put `def values: IndexedSeq[T] = findValues` here since ClpEnum is not a class!

    /** The string constructor method that enables us to create the case object instance by name. */
    def apply(name: String): T = withNameInsensitive(name)
  }

  /** Gets a description for use in a generating a plot.  Returns the `<sample> / <library>` if there is only one
    * sample and library, otherwise the input file name. */
  def plotDescription(reader: SamSource, input: FilePath): String = {
    val samples   = reader.readGroups.map(_.getSample).distinct
    val libraries = reader.readGroups.map(_.getLibrary).distinct
    (samples, libraries) match {
      case (Seq(sample), Seq(library)) => s"$sample / $library"
      case _                           => PathUtil.basename(input, trimExt=true).toString
    }
  }
}
