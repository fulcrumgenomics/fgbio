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

import scala.tools.nsc.Settings
import scala.tools.nsc.doc.base.comment._
import scala.tools.nsc.doc.html.Doclet
import scala.tools.nsc.doc.model.DocTemplateEntity
import scala.tools.nsc.reporters.ConsoleReporter

/** Case class to capture information about a field/column in a metrics class/file. */
case class ColumnDescription(name: String, typ: String, description: String)

/** Case class to capture information about a metrics class/file. */
case class MetricDescription(name: String, description: String, columns: Seq[ColumnDescription]) {
  def summary: String = description.takeWhile(_ != '.').replace('\n', ' ')
}

/**
  * Custom scaladoc Doclet for rendering the documentation for [[com.fulcrumgenomics.util.Metric]] classes into
  * MarkDown for display on the fgbio website.
  */
class FgMetricsDoclet extends Doclet(reporter = new ConsoleReporter(new Settings())) {
  /**
    * Main entry point for the doclet.  Scans for documentation for the metrics classes and
    * renders it into MarkDown.
    */
  override def generateImpl(): Unit = {
    // The MarkDown file to be written
    val md  = Paths.get(this.universe.settings.outdir.value).resolve("metrics.md")
    val out = new PrintStream(md.toFile)

    out.println(
      s"""
         |# fgbio Metrics Descriptions
         |
         |This page contains descriptions of all metrics produced by all fgbio tools.  Within the descriptions
         |the type of each field/column is given, including two commonly used types:
         |
         |* `Count` is an integer representing the count of some item
         |* `Proportion` is a real number with a value between 0 and 1 representing a proportion or fraction
         |
         |## Table of Contents
         |
         ||Metric Type|Description|
         ||-----------|-----------|""".stripMargin
    )

    metrics.foreach { m =>
      out.println(s"|[${m.name}](#${toLinkTarget(m.name)})|${m.summary}|")
    }

    out.println("\n## Metric File Descriptions")

    metrics.foreach { m =>
      out.println()
      out.println(s"\n### ${m.name}\n\n${m.description}\n")

      // The table of columns
      out.println("|Column|Type|Description|")
      out.println("|------|----|-----------|")
      m.columns.foreach { c =>
        out.println(s"|${c.name}|${c.typ}|${c.description}|")
      }
    }
  }

  /** Locates the metrics documentation templates and turns them into simple case classes with comments as markdown. */
  private lazy val metrics: Seq[MetricDescription] = {
    def simplify(name: String) = if (name.indexOf('.') > 0) name.substring(name.lastIndexOf('.') + 1) else name

    findMetricsClasses.map{ template =>
      val name        = template.name
      val description = template.comment.map(c => renderBody(c.body)).getOrElse("")
      val columns     = template.constructors.find(_.isPrimary) match {
        case None              => Seq.empty
        case Some(constructor) =>
          val comments = constructor.comment.map(c => c.valueParams).getOrElse(Map.empty[String,Body])
          constructor.valueParams.flatten.map { param =>
            val d    = comments.get(param.name).map(renderBody).getOrElse("").replace('\n', ' ')
            val desc = d.take(1).toUpperCase + d.drop(1)
            ColumnDescription(name=param.name, typ=simplify(param.resultType.name), description=desc)
          }
      }

      MetricDescription(name=name, description=description, columns=columns)
    }.sortBy(_.name)
  }

  /** Finds the [[scala.tools.nsc.doc.model.DocTemplateEntity]] instances that correspond to subclasses of [[Metric]] */
  private def findMetricsClasses: List[DocTemplateEntity] = {
    def find(template: DocTemplateEntity): List[DocTemplateEntity] = {
      template :: template.templates.collect { case d: DocTemplateEntity => find(d) }.flatten
    }

    find(universe.rootPackage)
      .filter(d => d.isClass && !d.isAbstract)
      .filter(d => d.parentTypes.exists { case (template, typ) => template.toString == classOf[Metric].getName })
  }

  /** Take the body of a scaladoc comment and renders it into MarkDown. */
  private def renderBody(body: Body): String = {
    val buffer = new StringBuilder

    // Takes a block element and renders it into MarkDown and writes it into the buffer
    def renderBlock(block: Block, indent: String): Unit = {
      (block: @unchecked) match {
        case para:  Paragraph      => render(para.text)
        case dlist: DefinitionList => () // TODO
        case hr:    HorizontalRule => () // TODO
        case olist: OrderedList    => () // TODO
        case title: Title          => buffer.append("#" * title.level).append(" "); render(title.text); buffer.append("\n\n")
        case ulist: UnorderedList  => () // TODO
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

  /** Turns the text from a heading into the text to use as a link target. */
  private def toLinkTarget(heading: String): String = heading.toLowerCase.replace(' ', '-')
}
