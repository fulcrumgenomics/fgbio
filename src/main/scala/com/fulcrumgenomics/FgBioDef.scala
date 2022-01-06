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

import scala.math.Ordering

/** FgBioDef that is just an extension of the commons CommonsDef. Here in case
  * there's a need for FgBio specific defs later.
  */
object FgBioDef extends CommonsDef {

  /** Extends this trait in your enumeration object to enable the case objects to be created on the command line.
    * You should implement the [[enumeratum.Enum#values]] method in your object, for example:
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


  /** Stores information about the coordinates of the alternate locus
    *
    * @param refName the name of the reference sequence (or chromosome)
    * @param start position of the genomic range (1-based)
    * @param end position of the gneomic range (1-based inclusive).
    */
  case class GenomicRange(refName: String, start: Int, end: Int) {
    override def toString: String = f"$refName:$start-$end"
  }

  object GenomicRange {
    /** Parses a genomic range of the form `<chr>:<start>-<end>` or `<chr>:<start>`, assuming the range is 1-based
      * inclusive (closed-ended).  If only a start position is included the range includes just the single base.
      * Coordinates may include commas in the numbers.
      *
      * @param range the string to parse.
      */
    def apply(range: String): GenomicRange = {
      val (refName: String, rest: String) = range.lastIndexOf(':') match {
        case -1  => throw new IllegalArgumentException(f"Missing colon in genomic range: '$range'")
        case idx =>
          val (left, right) = range.splitAt(idx)
          (left, right.drop(1).replace(",", ""))
      }

      val (start, end) = rest.indexOf('-') match {
        case -1  => (rest, rest)
        case idx =>
          val (left, right) = rest.splitAt(idx)
          (left, right.drop(1))
      }

      GenomicRange(refName=refName, start=start.toInt, end=end.toInt)
    }
  }

  /** Parses a genomic range of the form `<chr>:<start>-<end>` or `<chr>:<start>`, assuming the range is 1-based
    * inclusive (closed-ended).  If only a start position is included the range includes just the single base.
    * Coordinates may include commas in the numbers.
    *
    * @param range the string to parse.
    * @deprecated use GenomicRange(string) instead
    */
  @deprecated(message="Use GenomicRange() instead", since="1.3.0") def parseRange(range: String): GenomicRange = {
    GenomicRange(range)
  }

  // Developer note: move this to commons as scala 2.12 does not have this method on an [[IterableOnce]]/[[Iterator]]
  implicit class IterableOnceMinByOption[A](iter: Iterator[A]) {
    /** Finds the first element which yields the smallest value measured by function f.
      *
      * @param    cmp An ordering to be used for comparing elements.
      * @tparam   B The result type of the function f.
      * @param    f The measuring function.
      * @return an option value containing the first element of this $coll
      *         with the smallest value measured by function f
      *         with respect to the ordering `cmp`.
      */
    def minByOption[B](f: A => B)(implicit cmp: Ordering[B]): Option[A] = {
      if (iter.iterator.isEmpty) None
      else Some(iter.iterator.minBy(f)(cmp))
    }
  }
}
