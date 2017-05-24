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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.cmdline.FgBioTool
import com.fulcrumgenomics.sopt.Sopt.ClpMetadata
import com.fulcrumgenomics.sopt.{Sopt, arg, clp}
import com.fulcrumgenomics.util.Io

@clp(description="Generates the suite of per-tool MarkDown documents.")
class BuildToolDocs
( @arg(flag='o', doc="Output directory") val output: DirPath,
  @arg(flag='p', doc="The packages to document") val packages: Seq[String] = Seq("com.fulcrumgenomics"),
  @arg(flag='n', doc="The name of the tool chain") val name: String = "fgbio"
) extends InternalTool {
  
  override def execute(): Unit = {
    Io.mkdirs(output)
    val version = getClass.getPackage.getImplementationVersion
    val classes = Sopt.find[FgBioTool](packages=packages)
    logger.info(s"Found ${classes.size} tools to document.")
    
    val toolsByGroup = classes.map(c => Sopt.inspect(c)).groupBy(_.group.name)
    
    val indexHeader =
      s"""
        |---
        |title: $name tools
        |---
        |
        |# $name tools
        |
        |The following tools are available in $name version $version.
        |
      """.trim.stripMargin
    
    val builder = new StringBuilder
    builder.append(indexHeader)
    
    for (groupName <- toolsByGroup.keySet.toSeq.sorted) {
      val group = toolsByGroup(groupName).head.group
      
      builder.append(
        s"""
           |## ${group.name}
           |
           |${group.description}
           |
           ||Tool|Description|
           ||----|-----------|
         """.trim.stripMargin).append('\n')
  
      for (tool <- toolsByGroup(groupName).sortBy(_.name)) {
        val out = output.resolve(tool.name + ".md")
        Io.writeLines(out, Seq(documentTool(tool)))
        
        val shortDesc = tool.descriptionAsText.takeWhile(_ != '.').replace('\n', ' ')
        builder.append(s"|[${tool.name}](${tool.name}.md)|$shortDesc|\n")
      }
      
      builder.append("\n")
    }
    
    Io.writeLines(output.resolve("index.md"), Seq(builder.toString()))
  }
  
  def documentTool(tool: ClpMetadata): String = {
    val doc =
      s"""
         |---
         |title: ${tool.name}
         |---
         |
         |# ${tool.name}
         |
         |## Overview
         |**Group:** ${tool.group.name}
         |
         |${tool.description}
         |
         |## Arguments
         |
         ||Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
         ||----|----|----|-----------|---------|----------|----------------|
        """.stripMargin.trim
        
    val builder = new StringBuilder(doc).append('\n')
    tool.args.foreach { a => 
      val cols = Seq(
        a.name,
        a.flag.getOrElse(""),
        a.kind,
        a.description.lines.mkString(" "),
        if (a.minValues == 0) "Optional" else "Required",
        if (a.maxValues == Int.MaxValue) "Unlimited" else a.maxValues,
        a.defaultValues.mkString(", ")
      )
          
      cols.foreach(c => builder.append('|').append(c))
      builder.append('|').append('\n')
    }
    
    builder.toString()
  }
}
