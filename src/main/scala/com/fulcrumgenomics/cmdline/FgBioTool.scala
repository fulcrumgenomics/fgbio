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

import com.fulcrumgenomics.cmdline.FgBioMain.FailureException
import com.fulcrumgenomics.sopt.cmdline.ValidationException

/** All fgbio tools should extend this. */
// Todo: extend Lazy Logging?
trait FgBioTool {
  def execute(): Unit

  /** Fail with just an exit code. */
  def fail(exit: Int) = throw new FailureException(exit=exit)

  /** Fail with the default exit code and a message. */
  def fail(message: String) = throw new FailureException(message=Some(message))

  /** Fail with a specific error code and message. */
  def fail(exit: Int, message: String) = throw new FailureException(exit=exit, message=Some(message))

  /** Generates a new validation exception with the given message. */
  def invalid(message: String) = throw new ValidationException(message)

  /** Generates a validation exception if the test value is false. */
  def validate(test: Boolean, message: => String) = if (!test) throw new ValidationException(message)
}
