package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.UnitSpec

class ContigNameMappingTest extends UnitSpec {
  private val mappingLinesOneToOne = Seq(
    "NC_000001.10\tchr1",
    "NC_000002.11\tchr2",
    "NC_000003.11\tchr3",
    "NC_000004.11\tchr4"
  )
  private val mappingLinesOneToMany = Seq(
      "1\tchr1",
      "2\tchr2\tChr2",
      "3\tchr3",
      "3\tChr3",
  )

  private def mappingInput(mappingLines: Seq[String]) = {
      val path = makeTempFile("test.", ".txt")
      Io.writeLines(path, mappingLines)
      path
  }

  "ContigNameMappings" should "properly map 1:1 source target contigs" in {
      val srcToTarget = ContigNameMapping.parse(mappingInput(mappingLinesOneToOne))
      srcToTarget shouldBe Map("NC_000001.10" -> List("chr1"),
                               "NC_000002.11" -> List("chr2"),
                               "NC_000003.11" -> List("chr3"),
                               "NC_000004.11" -> List("chr4"))
  }

  it should "properly map one to many source target contigs" in {
      val srcToTarget = ContigNameMapping.parse(mappingInput(mappingLinesOneToMany))
      srcToTarget shouldBe Map("1" -> List("chr1"),
                               "2" -> List("chr2", "Chr2"),
                               "3" -> List("chr3", "Chr3"))
  }
}