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


package com.fulcrumgenomics.internal

import java.io.PrintStream
import java.nio.file.Paths

import com.fulcrumgenomics.util.Metric

import scala.tools.nsc.doc.base.comment._
import scala.tools.nsc.doc.html.Doclet
import scala.tools.nsc.doc.model.DocTemplateEntity



/** Case class to capture information about a field/column in a metrics class/file. */
case class ColumnDescription(name: String, typ: String, description: String)

/** Case class to capture information about a metrics class/file. */
case class MetricDescription(name: String, description: String, columns: Seq[ColumnDescription]) {
  def summary: String = description.takeWhile(_ != '.').replace('\n', ' ')
}

/**
  * Custom scaladoc Doclet for rendering the documentation for [[Metric]] classes into
  * MarkDown for display on the fgbio website.
  */
trait FgBioDoclet extends Doclet {

  /** Take the body of a scaladoc comment and renders it into MarkDown. */
  protected def renderBody(body: Body): String = {
    val buffer = new StringBuilder

    // Takes a block element and renders it into MarkDown and writes it into the buffer
    def renderBlock(block: Block, indent: String): Unit = {
      block match {
        case para:  Paragraph      => render(para.text)
        case dlist: DefinitionList => Unit // TODO
        case hr:    HorizontalRule => Unit // TODO
        case olist: OrderedList    => Unit // TODO
        case title: Title          => buffer.append("#" * title.level).append(" "); render(title.text); buffer.append("\n\n")
        case ulist: UnorderedList  => Unit // TODO
      }
    }

    // Takes an inline element and renders it into MarkDown and writes it into the buffer
    def render(inline: Inline): Unit = inline match {
      case bold:    Bold        => buffer.append("**"); render(bold.text); buffer.append("**")
      case chain:   Chain       => chain.items.foreach(render)
      case link:    EntityLink  => render(link.title) // TODO: better handling of entity links?
      case tag:     HtmlTag     => buffer.append(tag.data)
      case italic:  Italic      => buffer.append("*"); render(italic.text); buffer.append("__")
      case link:    Link        => buffer.append("[").append(link.target).append("]("); render(link.title); buffer.append(")")
      case mono:    Monospace   => buffer.append("`"); render(mono.text); buffer.append("`")
      case sub:     Subscript   =>buffer.append("<sub>"); render(sub.text); buffer.append("</sub>")
      case summary: Summary     => render(summary.text)
      case supe:    Superscript => buffer.append("<sup>"); render(supe.text); buffer.append("</sup>")
      case text:    Text        => buffer.append(text.text)
      case under:   Underline   => buffer.append("__"); render(under.text); buffer.append("__")
    }

    body.blocks.foreach(renderBlock(_, ""))
    buffer.toString()
  }

  /** Simplifies the name by removing any characters prior to the last period. */
  protected def simplify(name: String) = if (name.indexOf('.') > 0) name.substring(name.lastIndexOf('.') + 1) else name

  /** Turns the text from a heading into the text to use as a link target. */
  protected def toLinkTarget(heading: String): String = heading.toLowerCase.replace(' ', '-')
}

