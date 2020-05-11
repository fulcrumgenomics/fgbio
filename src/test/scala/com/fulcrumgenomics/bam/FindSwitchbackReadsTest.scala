/*
 * The MIT License
 *
 * Copyright (c) 2018 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.bam

import com.fulcrumgenomics.bam.FindSwitchbackReads.SwitchMetric
import com.fulcrumgenomics.bam.api.SamOrder
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.fasta.Converters.FromSAMSequenceDictionary
import com.fulcrumgenomics.fasta.ReferenceSequenceIterator
import com.fulcrumgenomics.testing.SamBuilder.{Minus, Plus}
import com.fulcrumgenomics.testing.{ReferenceSetBuilder, SamBuilder, UnitSpec}
import com.fulcrumgenomics.util.Metric
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import org.scalatest.OptionValues

class FindSwitchbackReadsTest extends UnitSpec with OptionValues {
  private val (ref, refPath) = {
    val builder = new ReferenceSetBuilder(assembly=Some("TF"))
    val seq = builder.add("tf")
    // Following sequence is 70bp per line
    """
      |TGGCACCGAGCGAGCCGCGATGACAATGGCTGCATTGTGCTTCATGTCCCTTCCCATCAACATTTCTGTG
      |CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAAAATGAGGTGATCAGTGGGACGAGTAAGGA
      |AGGGGGGTTGGGAGAGGGGCGATTGGGCAACCCGGCTGCACAAACACGGGAGGTCAAAGATTGCGCCCAG
      |CCCGCCCAGGCCGGGAATGGAATAAAGGGACGCGGGGCGCCGGAGGCTGCACAGAAGCGAGTCCGACTGT
      |GCTCGCTGCTCAGCGCCGCACCCGGAAGATGAGGCTCGCCGTGGGAGCCCTGCTGGTCTGCGCCGTCCTG
      |GGGCTGTGTCTGGCTGTCCCTGATAAAACTGTGAGATGGTGTGCAGTGTCGGAGCATGAGGCCACTAAGT
      |GCCAGAGTTTCCGCGACCATATGAAAAGCGTCATTCCATCCGATGGTCCCAGTGTTGCTTGTGTGAAGAA
      |AGCCTCCTACCTTGATTGCATCAGGGCCATTGCGGCAAACGAAGCGGATGCTGTGACACTGGATGCAGGT
      |TTGGTGTATGATGCTTACCTGGCTCCCAATAACCTGAAGCCTGTGGTGGCAGAGTTCTATGGGTCAAAAG
      |AGGATCCACAGACTTTCTATTATGCTGTTGCTGTGGTGAAGAAGGATAGTGGCTTCCAGATGAACCAGCT
      |TCGAGGCAAGAAGTCCTGCCACACGGGTCTAGGCAGGTCCGCTGGGTGGAACATCCCCATAGGCTTACTT
      |TACTGTGACTTACCTGAGCCACGTAAACCTCTTGAGAAAGCAGTGGCCAATTTCTTCTCGGGCAGCTGTG
      |CCCCTTGTGCGGATGGGACGGACTTCCCCCAGCTGTGTCAACTGTGTCCAGGGTGTGGCTGCTCCACCCT
      |TAACCAATACTTCGGCTACTCGGGAGCCTTCAAGTGTCTGAAGGATGGTGCTGGGGATGTGGCCTTTGTC
      |AAGCACTCGACTATATTTGAGAACTTGGCAAACAAGGCTGACAGGGACCAGTATGAGCTGCTTTGCCTGG
      |ACAACACCCGGAAGCCGGTAGATGAATACAAGGACTGCCACTTGGCCCAGGTCCCTTCTCATACCGTCGT
      |GGCCCGAAGTATGGGCGGCAAGGAGGACTTGATCTGGGAGCTTCTCAACCAGGCCCAGGAACATTTTGGC
      |AAAGACAAATCAAAAGAATTCCAACTATTCAGCTCTCCTCATGGGAAGGACCTGCTGTTTAAGGACTCTG
      |CCCACGGGTTTTTAAAAGTCCCCCCCAGGATGGATGCCAAGATGTACCTGGGCTATGAGTATGTCACTGC
      |CATCCGGAATCTACGGGAAGGCACATGCCCAGAAGCCCCAACAGATGAATGCAAGCCTGTGAAGTGGTGT
      |GCGCTGAGCCACCACGAGAGGCTCAAGTGTGATGAGTGGAGTGTTAACAGTGTAGGGAAAATAGAGTGTG
      |TATCAGCAGAGACCACCGAAGACTGCATCGCCAAGATCATGAATGGAGAAGCTGATGCCATGAGCTTGGA
      |TGGAGGGTTTGTCTACATAGCGGGCAAGTGTGGTCTGGTGCCTGTCTTGGCAGAAAACTACAATAAGAGC
      |GATAATTGTGAGGATACACCAGAGGCAGGGTATTTTGCTGTAGCAGTGGTGAAGAAATCAGCTTCTGACC
      |TCACCTGGGACAATCTGAAAGGCAAGAAGTCCTGCCATACGGCAGTTGGCAGAACCGCTGGCTGGAACAT
      |CCCCATGGGCCTGCTCTACAATAAGATCAACCACTGCAGATTTGATGAATTTTTCAGTGAAGGTTGTGCC
      |CCTGGGTCTAAGAAAGACTCCAGTCTCTGTAAGCTGTGTATGGGCTCAGGCCTAAACCTGTGTGAACCCA
      |ACAACAAAGAGGGATACTACGGCTACACAGGCGCTTTCAGGTGTCTGGTTGAGAAGGGAGATGTGGCCTT
      |TGTGAAACACCAGACTGTCCCACAGAACACTGGGGGAAAAAACCCTGATCCATGGGCTAAGAATCTGAAT
      |GAAAAAGACTATGAGTTGCTGTGCCTTGATGGTACCAGGAAACCTGTGGAGGAGTATGCGAACTGCCACC
      |TGGCCAGAGCCCCGAATCACGCTGTGGTCACACGGAAAGATAAGGAAGCTTGCGTCCACAAGATATTACG
      |TCAACAGCAGCACCTATTTGGAAGCAACGTAACTGACTGCTCGGGCAACTTTTGTTTGTTCCGGTCGGAA
      |ACCAAGGACCTTCTGTTCAGAGATGACACAGTATGTTTGGCCAAACTTCATGACAGAAACACATATGAAA
      |AATACTTAGGAGAAGAATATGTCAAGGCTGTTGGTAACCTGAGAAAATGCTCCACCTCATCACTCCTGGA
      |AGCCTGCACTTTCCGTAGACCTTAAAATCTCAGAGGTAGGGCTGCCACCAAGGTGAAGATGGGAACGCAG
      |ATGATCCATGAGTTTGCCCTGGTTTCACTGGCCCAAGTGGTTTGTGCTAACCACGTCTGTCTTCACAGCT
      |CTGTGTTGCCATGTGTGCTGAACAAAAAATAAAAATTATTATTGATTTTATATTTCAAAAACTCCATTCT
      |TTCCTAAATATTTTCAACAAAGGATTTCTTTATGCATTCTGCCTAAATACCTATGCAACTGAGCCCTTCC
      |TTCTCAGCTCAAGATTCGTCTGGTCTTTCCCTACAGCTTTGTGTGTGCCATGGCCACATCTCCTGGGTAC
      |AGTTCAAGGAGACATCTTTTCTAAAAGGGTCTGCGTGATCATTAAAATATAATCAAATGTAAAAAAAAAA
      |AAAAAAAA
    """.stripMargin.linesIterator.foreach(seq.add(_))

    builder.add("chr2").add("ACGTACGT", 100)

    val refPath = builder.toTempFile()
    val ref = ReferenceSequenceFileFactory.getReferenceSequenceFile(refPath)
    (ref, refPath)
  }

  private val refMap = ReferenceSequenceIterator(refPath, stripComments=true).map(r => r.getName -> r).toMap

  "FindSwitchbackReads.findReadBasedSwitchback" should "not find switchbacks when there is too little clipping" in {
    val builder = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))

    val noClip = Template(builder.addFrag(contig=0, start=71, cigar="50M", bases="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAAAATGAGGTG").iterator)
    val noClipHit = FindSwitchbackReads.findReadBasedSwitchback(refMap, noClip, minLength=5, maxOffset=20, maxErrorRate=0.1)
    noClipHit shouldBe None

    val fourClip = Template(builder.addFrag(contig=0, start=71, cigar="4S46M", bases="ccagCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAAAATGA").iterator)
    val fourClipHit = FindSwitchbackReads.findReadBasedSwitchback(refMap, noClip, minLength=5, maxOffset=20, maxErrorRate=0.1)
    fourClipHit shouldBe None
  }

  it should "find a switchback at the start of a plus strand read" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val template = Template(builder.addFrag(contig=0, start=71, cigar="10S40M", bases="aggagtccagCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGA").iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit.value.length shouldBe 10
    hit.value.offset shouldBe 0
  }

  it should "not find a switchback at the start of a plus strand read when the error rate is too high" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val template = Template(builder.addFrag(contig=0, start=71, cigar="10S40M", bases="aCCagtccagCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGA").iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit shouldBe None
  }

  it should "find a switchback at the start of a plus strand read with a large offset" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val template = Template(builder.addFrag(contig=0, start=71, cigar="10S40M", bases="ttgatgggaaCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGA").iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit.value.length shouldBe 10
    hit.value.offset shouldBe 20
  }

  it should "not find a switchback at the start of a plus strand read with an offset that is too large (25bp>20bp)" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val template = Template(builder.addFrag(contig=0, start=71, cigar="10S40M", bases="gggaagggacCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGA").iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit shouldBe None
  }

  it should "not find a switchback in a plus strand read with soft-clipping at the _end_" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val template = Template(builder.addFrag(contig=0, start=71, cigar="40M10S", bases="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAtctgagctct").iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit shouldBe None
  }

  it should "find a switchback in a negative strand read with soft-clipping at the end" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val template = Template(builder.addFrag(contig=0, start=71, cigar="40M10S", strand=Minus, bases="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAtctgagctct").iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit.value.length shouldBe 10
    hit.value.offset shouldBe 0
  }

  it should "calculate offset correctly for negative strand reads" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val plus10  = Template(builder.addFrag(contig=0, start=71, cigar="40M10S", strand=Minus, bases="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAcacctcattt").iterator)
    val minus10 = Template(builder.addFrag(contig=0, start=71, cigar="40M10S", strand=Minus, bases="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAggagacgacc").iterator)
    val plusHit  = FindSwitchbackReads.findReadBasedSwitchback(refMap, plus10, minLength=5, maxOffset=20, maxErrorRate=0.1)
    val minusHit = FindSwitchbackReads.findReadBasedSwitchback(refMap, minus10, minLength=5, maxOffset=20, maxErrorRate=0.1)
    plusHit.value.offset shouldBe 10
    minusHit.value.offset shouldBe -10
  }

  it should "find a switchback with errors in the sequence at the start of a read" in {
    val builder  = new SamBuilder(readLength=100, sd=Some(ref.getSequenceDictionary.fromSam))
    val rec = builder.addFrag(contig=0, start=71, cigar="50S50M", bases="cacctcattttctgagctctggagacgacccgcgagtggaaggagtccagCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAAAATGAGGTG").get
    val template = Template(Iterator(rec))

    // No Errors
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe true

    rec.bases(0)  = 'g'
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe true
    rec.bases(5)  = 't'
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe true
    rec.bases(10) = 'a'
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe true
    rec.bases(15) = 'c'
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe true
    rec.bases(20) = 't'
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe true
    rec.bases(25) = 'g'
    FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1).isDefined shouldBe false
  }

  it should "find a switchback in a read when given a template with paired reads" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val recs     = builder.addPair(contig=0, start1=1,  cigar1="50M",    strand1=Plus,  bases1="TGGCACCGAGCGAGCCGCGATGACAATGGCTGCATTGTGCTTCATGTCCC",
                                             start2=71, cigar2="40M10S", strand2=Minus, bases2="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAtctgagctct")
    val template = Template(recs.iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit.value.length shouldBe 10
    hit.value.offset shouldBe 0
    hit.value.read shouldBe recs.find(_.secondOfPair)
  }

  it should "not find switchbacks in secondary or supplementary records" in {
    val builder  = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    builder.addFrag(name="q1", contig=0, start=71, cigar="50M", strand=Plus, bases="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAAAATGAGGTG")
    builder.addFrag(name="q1", contig=0, start=71, cigar="10S40M", strand=Plus, bases="aggagtccagCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGA").foreach(_.secondary = true)
    builder.addFrag(name="q1", contig=0, start=71, cigar="12S38M", strand=Plus, bases="gaaggagtccagCTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCA").foreach(_.supplementary = true)
    val template = Template(builder.iterator)
    val hit      = FindSwitchbackReads.findReadBasedSwitchback(refMap, template, minLength=5, maxOffset=20, maxErrorRate=0.1)
    hit shouldBe None
  }

  "FindSwitchbackReads.findTandemSwitchback" should "call FF and RR templates with small enough gaps as switchbacks" in {
    val builder = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val ff = Template(builder.addPair(contig=0, start1=100, start2=160, strand1=Plus, strand2=Plus).iterator)
    val rr = Template(builder.addPair(contig=0, start1=100, start2=160, strand1=Minus, strand2=Minus).iterator)

    FindSwitchbackReads.findTandemSwitchback(ff, maxGap=200).value.gap shouldBe 10
    FindSwitchbackReads.findTandemSwitchback(rr, maxGap=200).value.gap shouldBe 10
  }

  it should "not call tandem pairs switchbacks if the gap is too big" in {
    val readLen = 250
    val start   = 500
    val builder = new SamBuilder(readLength=readLen, sd=Some(ref.getSequenceDictionary.fromSam))

    (Seq(-201, 201) ++ Range(-200, 200, step=20)).foreach { gap =>
      val ff  = Template(builder.addPair(contig=0, start1=start, start2=start+readLen+gap, strand1=Plus, strand2=Plus).iterator)
      val hit = FindSwitchbackReads.findTandemSwitchback(ff, maxGap=200)

      if (gap <= 200 && gap >= -200) hit.value.gap shouldBe gap
      else hit shouldBe None
    }
  }

  it should "not call tandem pairs switchbacks if the reads are on different chromosomes" in {
    val builder = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val ff = Template(builder.addPair(contig=0, start1=100, start2=160, strand1=Plus, strand2=Plus).iterator)
    FindSwitchbackReads.findTandemSwitchback(ff, maxGap=200) shouldBe defined

    ff.r2.value.refIndex = 1
    ff.r1.value.mateRefIndex = 1
    FindSwitchbackReads.findTandemSwitchback(ff, maxGap=200) shouldBe None
  }

  it should "not call FR or RF pairs switchbacks" in {
    val builder = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam))
    val ff = Template(builder.addPair(contig=0, start1=100, start2=160, strand1=Plus, strand2=Plus).iterator)
    val fr = Template(builder.addPair(contig=0, start1=100, start2=160, strand1=Plus, strand2=Minus).iterator)
    val rf = Template(builder.addPair(contig=0, start1=100, start2=160, strand1=Minus, strand2=Plus).iterator)

    FindSwitchbackReads.findTandemSwitchback(ff, maxGap=1000) shouldBe defined
    FindSwitchbackReads.findTandemSwitchback(fr, maxGap=1000) shouldBe None
    FindSwitchbackReads.findTandemSwitchback(rf, maxGap=1000) shouldBe None
  }

  "FindSwitchbackReads" should "run end to end and produce a tagged BAM and metrics" in {
     val builder = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam), sort=Some(SamOrder.Queryname))
    builder.rg.setSample("switchback_sample")
    builder.rg.setLibrary("switchback_library")

    // Add 20 "normal" read pairs
    Range.inclusive(10, 200, step=10).foreach(s => builder.addPair(name=s"fr$s", start1=s, start2=s+200))

    // Add an unmapped pair
    builder.addPair(name="u1", unmapped1=true, unmapped2=true)

    // Add 6 FF pairs
    Range.inclusive(100, 600, step=100).foreach(s => builder.addPair(name=s"ff$s", start1=s, start2=s+60, strand1=Plus, strand2=Plus))

    // Add 5 RR pairs
    Range.inclusive(1100, 1500, step=100).foreach(s => builder.addPair(name=s"rr$s", start1=s, start2=s+80, strand1=Minus, strand2=Minus))

    // Add a switchback with length=10, offset=0 in R2, with a pair of supplementary alignments
    builder.addPair(name="q1", contig=0, start1=1,  cigar1="50M",    strand1=Plus,  bases1="TGGCACCGAGCGAGCCGCGATGACAATGGCTGCATTGTGCTTCATGTCCC",
                                         start2=71, cigar2="40M10S", strand2=Minus, bases2="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAtctgagctct")
    builder.addPair(name="q1", contig=0, start1=580, start2=670, bases1="TGGCACCGAGCGAGCCGCGATGACAATGGCTGCATTGTGCTTCATGTCCC",
      bases2="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAtctgagctct").foreach(_.supplementary = true)

    // Add a switchback with length=20, offset=20 in R1 with a pair of secondary alignments
    builder.addPair(name="q2", contig=0, start1=91,  cigar1="20S30M", strand1=Plus,  bases1="cgcgagtggaaggagtccagGGTCGTCTCCAGAGCTCAGAAAATGAGGTG",
                                         start2=211, cigar2="50M",    strand2=Minus, bases2="CCCGCCCAGGCCGGGAATGGAATAAAGGGACGCGGGGCGCCGGAGGCTGC")
    builder.addPair(name="q2", start1=754, start2=153, bases1="cgcgagtggaaggagtccagGGTCGTCTCCAGAGCTCAGAAAATGAGGTG",
      bases2="CCCGCCCAGGCCGGGAATGGAATAAAGGGACGCGGGGCGCCGGAGGCTGC").foreach(_.secondary = true)

    val in = builder.toTempFile()
    val out = makeTempFile("find_switchbacks.", ".bam")
    val metricBase = makeTempFile("find_switcbacks.", ".metrics")
    new FindSwitchbackReads(input=in, output=out, metrics=Some(metricBase), ref=refPath).execute()

    // Assert some things about the BAM file
    val recs = readBamRecs(out)
    recs.foreach { rec =>

      rec.name match {
        case "u1" =>
          // Unmapped reads should be untouched
          rec.unmapped shouldBe true
          rec.get[String]("OA").isDefined shouldBe false
          rec.get[String](FindSwitchbackReads.SwitchTag).isDefined shouldBe false
        case n if n.startsWith("fr") =>
          // Regular FR reads with no swithbacks should be untouched
          rec.mapped shouldBe true
          rec.get[String]("OA").isDefined shouldBe false
          rec.get[String](FindSwitchbackReads.SwitchTag).isDefined shouldBe false
        case n if n.startsWith("q") =>
          // Reads with read-based switchbacks should be unmapped and secondary+supplementary reads removed
          rec.secondary shouldBe false
          rec.supplementary shouldBe false
          rec.mapped shouldBe false
          rec.get[String]("OA").value should not be empty
          rec.get[String](FindSwitchbackReads.SwitchTag).value.startsWith("r,") shouldBe true
        case n if n.startsWith("ff") || n.startsWith("rr") =>
          // Tandem reads in the size range should be unmapped
          rec.mapped shouldBe false
          rec.get[String]("OA").value should not be empty
          rec.get[String](FindSwitchbackReads.SwitchTag).value.startsWith("t,") shouldBe true
      }
    }

    val ms = Metric.read[SwitchMetric](PathUtil.pathTo(s"${metricBase}.summary.txt")).head
    ms.sample                   shouldBe "switchback_sample"
    ms.library                  shouldBe "switchback_library"
    ms.templates                shouldBe 20+1+6+5+2
    ms.aligned_templates        shouldBe 20+6+5+2
    ms.switchback_templates     shouldBe 6+5+2
    ms.fraction_switchbacks     shouldBe (6+5+2)/(20+6+5+2).toDouble +- 0.001
    ms.read_based_switchbacks   shouldBe 2
    ms.mean_length              shouldBe 15.0
    ms.mean_offset              shouldBe 10.0
    ms.tandem_based_switchbacks shouldBe 11
    ms.mean_gap                 shouldBe (6*10 + 5*30) / 11.toDouble +- 0.001
  }

  it should "run end to end and not unmap reads when dont-unmap is set" in {
    val builder = new SamBuilder(readLength=50, sd=Some(ref.getSequenceDictionary.fromSam), sort=Some(SamOrder.Queryname))

    // Add an FF pair
    builder.addPair(name=s"ff", start1=500, start2=560, strand1=Plus, strand2=Plus)

    // Add a switchback with length=10, offset=0 in R2
    builder.addPair(name="q1", contig=0, start1=1,  cigar1="50M",    strand1=Plus,  bases1="TGGCACCGAGCGAGCCGCGATGACAATGGCTGCATTGTGCTTCATGTCCC",
                                         start2=71, cigar2="40M10S", strand2=Minus, bases2="CTGGACTCCTTCCACTCGCGGGTCGTCTCCAGAGCTCAGAtctgagctct")

    val in = builder.toTempFile()
    val out = makeTempFile("find_switchbacks.", ".bam")
    new FindSwitchbackReads(input=in, output=out, metrics=None, ref=refPath, dontUnmap=true).execute()
    val recs = readBamRecs(out)
    recs.foreach { rec =>
      rec.mapped shouldBe true
      rec.get[String]("OA").isDefined shouldBe false
      rec.get[String](FindSwitchbackReads.SwitchTag).isDefined shouldBe true
    }
  }
}
