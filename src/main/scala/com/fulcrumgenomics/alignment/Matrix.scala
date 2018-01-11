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

package com.fulcrumgenomics.alignment

import com.fulcrumgenomics.FgBioDef._

import scala.reflect.ClassTag

/** Factory method(s) for Matrices. */
object Matrix {
  /**
    * Creates a matrix of the desired size
    * @param rows the number of rows in the matrix
    * @param columns the number of columns in the matrix
    * @tparam A the type of element the matrix will hold
    * @return the newly created Matrix
    */
  def apply[A: ClassTag](rows: Int, columns: Int): Matrix[A] = {
    require(rows >= 0, s"Requested rows ($rows) is less than 0.")
    require(columns >= 0, s"Requested columns ($columns) is less than 0.")
    new LinearMatrix(rows, columns)
  }
}


/**
  * Defines methods applicable to a 2D matrix.  Matrix cells are accessed with a pair of
  * coordinates, the first of which is the row offset starting at 0, and the second of
  * which is the column offset, starting at 0.
  * */
trait Matrix[A] {
  /** Access a cell in the matrix. */
  def apply(i: Int, j: Int): A

  /** Updates a cell in the matrix. */
  def update(i: Int, j: Int, value: A): Unit

  /**
    * Returns the dimensions of the matrix as a tuple of [[Int]]s. E.g. for a 2x2 matrix this method will
    * return (2,2).  Note that, like an array, since the coordinates into the matrix are zero-based,
    * the bottom-right cell in a 2x2 matrix is accessed using matrix(1,1).
    */
  def dimensions: (Int, Int)

  /** Returns the total number of cells in the matrix. */
  def size: Int = {
    val (x,y) = dimensions
    x * y
  }

  /** Provides a simplistic text dump of the matrix, mostly for debugging. */
  override def toString: String = {
    val (x,y) = dimensions
    Range(0, x).map { i => Range(0, y).map { j => String.valueOf(this(i,j)) }.mkString("\t") }.mkString("\n")
  }
}


/** Implements a matrix using a single linear array. */
class LinearMatrix[A: ClassTag](val x: Int, val y: Int) extends Matrix[A] {
  private val data = new Array[A](x * y)

  /** Access a cell in the matrix. */
  override def apply(i: Int, j: Int): A = {
    if (i >= x || j >= y || i < 0 || j < 0) throw new NoSuchElementException(s"Illegal matrix indices: (${i}, ${j})")
    data(i*y + j)
  }

  /** Updates a cell in the matrix. */
  override def update(i: Int, j: Int, value: A): Unit = {
    if (i >= x || j >= y || i < 0 || j < 0) throw new NoSuchElementException(s"Illegal matrix indices: (${i}, ${j})")
    data(i*y + j) = value
  }

  /** Returns the dimensions of the matrix. */
  override def dimensions: (Int, Int) = (x, y)
}
