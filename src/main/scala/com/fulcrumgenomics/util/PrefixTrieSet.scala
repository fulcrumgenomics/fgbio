/*
 * The MIT License
 *
 * Copyright (c) 2025 Fulcrum Genomics
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
package com.fulcrumgenomics.util

class PrefixTrieSet(key: Option[Char] = None) {

  private var nodes  = new scala.collection.immutable.TreeMap[Char, PrefixTrieSet]
  private var numSeen = 0L

  def add(k: String): Long = {
    if (k.isEmpty) {
      this.numSeen += 1
      this.numSeen
    }
    else {
      val key = k.charAt(0)
      val node = this.nodes.get(key) match {
        case Some(node) => node
        case None       => {
          val node = new PrefixTrieSet(Some(key))
          this.nodes = this.nodes.updated(key, node)
          node
        }
      }
      node.add(k.substring(1))
    }
  }

  def contains(k: String): Boolean = {
    nodeFor(k).isDefined
  }

  def prefixes(): Iterator[String] = {
    val thisIter = if (this.numSeen > 0) Iterator("") else Iterator.empty
    thisIter ++ this.nodes.iterator.flatMap { case (char, node) =>
      node.prefixes().map { suffix => f"${char}${suffix}" }
    }
  }

  def values(): Iterator[(String, Long)] = {
    val thisIter = if (this.numSeen > 0) Iterator(("", this.numSeen)) else Iterator.empty
    thisIter ++ this.nodes.iterator.flatMap { case (char, node) =>
      node.values().map { case (suffix, count) =>
        (f"${char}${suffix}", count)
      }
    }
  }

  private def nodeFor(k: String): Option[PrefixTrieSet] =  {
    if (k.isEmpty) {
      if (this.numSeen > 0) Some(this)
      else None
    }
    else {
      this.nodes.get(k.charAt(0)).flatMap { node =>
        node.nodeFor(k.substring(1))
      }
    }
  }
}
