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

import com.fulcrumgenomics.commons.reflect.ReflectionUtil
import com.fulcrumgenomics.umi.ConsensusTags
import com.fulcrumgenomics.umi.ConsensusTags.PerBase
import com.fulcrumgenomics.umi.ConsensusTags.PerRead
import com.fulcrumgenomics.util.Metric

import scala.tools.nsc.doc.base.comment._
import scala.tools.nsc.doc.html.Doclet
import scala.tools.nsc.doc.model.DocTemplateEntity

/**
  * Custom scaladoc Doclet for rendering the documentation for SAM tags produced by consensus calling tools described in
  * [[ConsensusTags]] into MarkDown for display on the fgbio website.
  */
class FgConsensusTagsDoclet extends FgBioDoclet {
  /**
    * Main entry point for the doclet.  Scans for documentation for the tags and
    * renders it into MarkDown.
    */
  override def generateImpl(): Unit = {
    // The MarkDown file to be written
    val md  = Paths.get(this.universe.settings.outdir.value).resolve("consensustags.md")
    val out = new PrintStream(md.toFile)

    out.println(
      s"""
         |# fgbio Consensus Tag Descriptions
         |
         |This page contains descriptions of all SAM tags produced by fgbio tools that call consensus reads:
         |- [CallDuplexConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html)
         |- [CallMolecularConsensusReads](http://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html)
         |
         |Standard SAM tags can be found in the [hts-spec SAM tags document](https://samtools.github.io/hts-specs/SAMtags.pdf).
         |
         |Within the descriptions the type of each field/column is given.
         |
         |## Table of Contents
         |
         ||Tag Group|Description|
         ||-----------|-----------|""".stripMargin
    )

    metrics.foreach { m =>
      out.println(s"|[${m.name}](#${toLinkTarget(m.name)})|${m.summary}|")
    }

    out.println("\n## SAM Tag Descriptions")

    metrics.foreach { m =>
      out.println()
      out.println(s"\n### ${m.name}\n\n${m.description}\n")

      // The table of columns
      out.println("|Name|Type|Description|")
      out.println("|------|----|-----------|")
      m.columns.foreach { c =>
        out.println(s"|${c.name}|${c.typ}|${c.description}|")
      }
    }
  }

  /** Converts camel case to a pretty format */
  private def formatColumnName(text: String): String = {
    text.map { c: Char => if (c.isUpper) " " + c else c.toString }
      .dropWhile(_ == ' ')
      .mkString
  }

  /** Locates the metrics documentation templates and turns them into simple case classes with comments as markdown. */
  private lazy val metrics: Seq[MetricDescription] = {

    def find(template: DocTemplateEntity): List[DocTemplateEntity] = {
      template :: template.templates.collect { case d: DocTemplateEntity => find(d) }.flatten
    }

    def findClass(obj: Any): DocTemplateEntity = {
      val name = obj.getClass.getName.replace('$', '.').replaceAll("""\.$""", "") // for objects
      val docs = find(universe.rootPackage)
        .filter(d => d.isObject)
        .filter(d => d.qualifiedName == name)
      docs match {
        case doc :: Nil => doc
        case _ => throw new IllegalArgumentException(s"Expected a single document for '${obj.getClass.getName}' but found ${docs.length}")
      }
    }

    val perBaseTags   = findClass(ConsensusTags.PerBase)
    val perReadTags   = findClass(ConsensusTags.PerRead)
    val consensusTags = findClass(ConsensusTags)

    Seq(consensusTags, perBaseTags, perReadTags).map { tags =>
      val name        = tags.name
      val description = tags.comment.map(c => renderBody(c.body)).getOrElse("")

      val columns = tags.values
        .filterNot(v => simplify(v.resultType.name).contains("Seq")) // for ConsensusTags.PerBase.AllPerBaseTags
        .map { v =>
          val d    = v.comment.map(c => renderBody(c.body)).getOrElse("").replace('\n', ' ')
          val desc = d.take(1).toUpperCase + d.drop(1)
          ColumnDescription(name=formatColumnName(v.name), typ=simplify(v.resultType.name), description=desc)
        }
      MetricDescription(name=formatColumnName(name), description=description, columns=columns)
    }
  }
}
