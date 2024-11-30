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

package com.fulcrumgenomics.vcf.validation

import com.fulcrumgenomics.commons.util.{LogLevel, Logger}
import com.fulcrumgenomics.vcf.api.{Genotype, Variant}


case class ValidationResult(message: String,
                            variant: Option[Variant] = None,
                            genotype: Option[Genotype] = None,
                            level: LogLevel = LogLevel.Error) {
  def fullMessage: String = {
    val builder = new StringBuilder()
    variant.foreach(v => builder.append(f"${v.chrom}:${v.pos} "))
    genotype.foreach(g => builder.append(f"${g.sample} "))
    builder.append(this.message)
    builder.toString()
  }

  def emit(print: String => Unit): Unit = print(this.message)
  def emit(logger: Logger): Unit = this._emit(logger)
  val _emit: Logger => Unit = this.level match {
    case LogLevel.Debug   => (logger: Logger) => logger.debug(this.fullMessage)
    case LogLevel.Info    => (logger: Logger) => logger.info(this.fullMessage)
    case LogLevel.Warning => (logger: Logger) => logger.warning(this.fullMessage)
    case LogLevel.Error   => (logger: Logger) => logger.error(this.fullMessage)
    case LogLevel.Fatal   => (logger: Logger) => logger.fatal(this.fullMessage)
  }
}

object ValidationResult {
  def debug(message: String, variant: Option[Variant] = None, genotype: Option[Genotype] = None): ValidationResult = {
    ValidationResult(variant=variant, message=message, genotype=genotype, level=LogLevel.Debug)
  }
  def debug(message: String, variant: Variant): ValidationResult = debug(message=message, variant=Some(variant))
  def debug(message: String, variant: Variant, genotype: Genotype): ValidationResult = debug(message=message, variant=Some(variant), genotype=Some(genotype))

  def info(message: String, variant: Option[Variant] = None, genotype: Option[Genotype] = None): ValidationResult = {
    ValidationResult(variant=variant, message=message, genotype=genotype, level=LogLevel.Info)
  }
  def info(message: String, variant: Variant): ValidationResult = info(message=message, variant=Some(variant))
  def info(message: String, variant: Variant, genotype: Genotype): ValidationResult = info(message=message, variant=Some(variant), genotype=Some(genotype))

  def warning(message: String, variant: Option[Variant] = None, genotype: Option[Genotype] = None): ValidationResult = {
    ValidationResult(variant=variant, message=message, genotype=genotype, level=LogLevel.Warning)
  }
  def warning(message: String, variant: Variant): ValidationResult = warning(message=message, variant=Some(variant))
  def warning(message: String, variant: Variant, genotype: Genotype): ValidationResult = warning(message=message, variant=Some(variant), genotype=Some(genotype))

  def error(message: String, variant: Option[Variant] = None, genotype: Option[Genotype] = None): ValidationResult = {
    ValidationResult(variant=variant, message=message, genotype=genotype, level=LogLevel.Error)
  }
  def error(message: String, variant: Variant): ValidationResult = error(message=message, variant=Some(variant))
  def error(message: String, variant: Variant, genotype: Genotype): ValidationResult = error(message=message, variant=Some(variant), genotype=Some(genotype))
}