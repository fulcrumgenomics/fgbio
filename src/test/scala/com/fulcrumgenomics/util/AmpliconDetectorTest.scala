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

package com.fulcrumgenomics.util

import com.fulcrumgenomics.testing.{SamBuilder, UnitSpec}
import htsjdk.samtools.SamPairUtil
import org.scalatest.OptionValues

class AmpliconDetectorTest extends UnitSpec with OptionValues {

  "AmpliconDetector.maxPrimerLength" should "return the maximum primer length" in {
    AmpliconDetector(
      amplicons = Seq(Amplicon("chr1", 1, 10, 20, 40)),
      slop      = 0,
      unclippedCoordinates = true
    ).maxPrimerLength shouldBe 21
    AmpliconDetector(
      amplicons = Seq(Amplicon("chr1", 1, 20, 30, 40)),
      slop      = 0,
      unclippedCoordinates = true
    ).maxPrimerLength shouldBe 20
    AmpliconDetector(
      amplicons = Seq(Amplicon("chr1", 1, 21, 20, 40)),
      slop      = 0,
      unclippedCoordinates = true
    ).maxPrimerLength shouldBe 21
  }

  "AmpliconDetector.findPrimer(refName: String, start: Int, end: Int, positiveStrand: Boolean)" should "not find any matches" in {
    val amplicon = Amplicon("chr1", 1, 20, 80, 100)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=0, unclippedCoordinates=true)
    // different refName
    detector.findPrimer(refName="chr2", start=1, end=20, positiveStrand=true).isEmpty shouldBe true
    // overlaps right primer, but positive strand
    detector.findPrimer(refName="chr1", start=80, end=100, positiveStrand=true).isEmpty shouldBe true
    // overlaps left primer, but negative strand
    detector.findPrimer(refName="chr1", start=1, end=20, positiveStrand=false).isEmpty shouldBe true
    // no overlaps
    detector.findPrimer(refName="chr1", start=101, end=120, positiveStrand=true).isEmpty shouldBe true
    detector.findPrimer(refName="chr1", start=101, end=120, positiveStrand=false).isEmpty shouldBe true
    // overlaps the amplicon, but not the primer
    detector.findPrimer(refName="chr1", start=21, end=79, positiveStrand=true).isEmpty shouldBe true
    detector.findPrimer(refName="chr1", start=21, end=79, positiveStrand=false).isEmpty shouldBe true
  }

  it should "find a match for overlaps" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=10, unclippedCoordinates=true)

    // left primer, exactly
    detector.findPrimer(refName="chr1", start=100, end=120, positiveStrand=true).value shouldBe amplicon

    // left primer, within slop
    detector.findPrimer(refName="chr1", start=90, end=120, positiveStrand=true).value shouldBe amplicon // exactly slop
    detector.findPrimer(refName="chr1", start=110, end=120, positiveStrand=true).value shouldBe amplicon // exactly slop
    detector.findPrimer(refName="chr1", start=100, end=100, positiveStrand=true).value shouldBe amplicon // end coordinate not considered
    detector.findPrimer(refName="chr1", start=89, end=120, positiveStrand=true).isEmpty shouldBe true // to prove slop boundary
    detector.findPrimer(refName="chr1", start=111, end=120, positiveStrand=true).isEmpty shouldBe true // to prove slop boundary

    // right primer, exactly
    detector.findPrimer(refName="chr1", start=180, end=200, positiveStrand=false).value shouldBe amplicon

    // right primer, within slop
    detector.findPrimer(refName="chr1", start=180, end=210, positiveStrand=false).value shouldBe amplicon // exactly slop
    detector.findPrimer(refName="chr1", start=180, end=190, positiveStrand=false).value shouldBe amplicon // exactly slop
    detector.findPrimer(refName="chr1", start=200, end=200, positiveStrand=false).value shouldBe amplicon // start coordinate not considered
    detector.findPrimer(refName="chr1", start=180, end=211, positiveStrand=false).isEmpty shouldBe true // to prove slop boundary
    detector.findPrimer(refName="chr1", start=180, end=189, positiveStrand=false).isEmpty shouldBe true // to prove slop boundary
  }

  it should "only return the hit with the smallest distance for multiple overlaps" in {
    val amplicons = IndexedSeq(Amplicon("chr1", 100, 120, 180, 200), Amplicon("chr1", 110, 130, 190, 210))
    val detector  = AmpliconDetector(amplicons=amplicons, slop=10, unclippedCoordinates=true)

    // amplicon #1 left primer, exactly
    detector.findPrimer(refName="chr1", start=100, end=120, positiveStrand=true).value shouldBe amplicons.head
    // amplicon #2 left primer, exactly
    detector.findPrimer(refName="chr1", start=110, end=130, positiveStrand=true).value shouldBe amplicons.last
    // closer to amplicon #1
    detector.findPrimer(refName="chr1", start=104, end=130, positiveStrand=true).value shouldBe amplicons.head
    // tie
    detector.findPrimer(refName="chr1", start=105, end=130, positiveStrand=true).nonEmpty shouldBe true
    // closer to amplicon #2
    detector.findPrimer(refName="chr1", start=106, end=130, positiveStrand=true).value shouldBe amplicons.last

    // amplicon #1 right primer, exactly
    detector.findPrimer(refName="chr1", start=180, end=200, positiveStrand=false).value shouldBe amplicons.head
    // amplicon #2 right primer, exactly
    detector.findPrimer(refName="chr1", start=190, end=210, positiveStrand=false).value shouldBe amplicons.last
    // closer to amplicon #1
    detector.findPrimer(refName="chr1", start=180, end=204, positiveStrand=false).value shouldBe amplicons.head
    // tie
    detector.findPrimer(refName="chr1", start=180, end=205, positiveStrand=false).nonEmpty shouldBe true
    // closer to amplicon #2
    detector.findPrimer(refName="chr1", start=180, end=206, positiveStrand=false).value shouldBe amplicons.last
  }

  "AmpliconDetector.findPrimer(rec: SamRecord)" should "find matches using the aligned/clipped coordinates" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=5, unclippedCoordinates=false)
    val builder  = new SamBuilder()

    // positive strand
    detector.findPrimer(rec=builder.addFrag(start=100, cigar="10S90M40H").value).value shouldBe amplicon
    detector.findPrimer(rec=builder.addFrag(start=95, cigar="10S90M40H").value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=94, cigar="10S90M40H").value).isEmpty shouldBe true // outside slop
    detector.findPrimer(rec=builder.addFrag(start=105, cigar="10S90M40H").value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=1056, cigar="10S90M40H").value).isEmpty shouldBe true // outside slop

    // negative strand
    detector.findPrimer(rec=builder.addFrag(start=110, cigar="10S90M40H", strand=SamBuilder.Minus).value).value shouldBe amplicon
    detector.findPrimer(rec=builder.addFrag(start=106, cigar="10S90M40H", strand=SamBuilder.Minus).value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=105, cigar="10S90M40H", strand=SamBuilder.Minus).value).isEmpty shouldBe true // outside slop
    detector.findPrimer(rec=builder.addFrag(start=116, cigar="10S90M40H", strand=SamBuilder.Minus).value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=117, cigar="10S90M40H", strand=SamBuilder.Minus).value).isEmpty shouldBe true // v slop
  }

  it should "find matches using the unclipped coordinates" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=5, unclippedCoordinates=true)
    val builder  = new SamBuilder()

    // positive strand
    detector.findPrimer(rec=builder.addFrag(start=110, cigar="10S90M40H").value).value shouldBe amplicon
    detector.findPrimer(rec=builder.addFrag(start=105, cigar="10S90M40H").value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=104, cigar="10S90M40H").value).isEmpty shouldBe true // outside slop
    detector.findPrimer(rec=builder.addFrag(start=115, cigar="10S90M40H").value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=116, cigar="10S90M40H").value).isEmpty shouldBe true // outside slop

    // negative strand
    detector.findPrimer(rec=builder.addFrag(start=71, cigar="10S90M40H", strand=SamBuilder.Minus).value).value shouldBe amplicon
    detector.findPrimer(rec=builder.addFrag(start=66, cigar="10S90M40H", strand=SamBuilder.Minus).value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=65, cigar="10S90M40H", strand=SamBuilder.Minus).value).isEmpty shouldBe true // outside slop
    detector.findPrimer(rec=builder.addFrag(start=76, cigar="10S90M40H", strand=SamBuilder.Minus).value).value shouldBe amplicon // within slop
    detector.findPrimer(rec=builder.addFrag(start=77, cigar="10S90M40H", strand=SamBuilder.Minus).value).isEmpty shouldBe true // v slop
  }

  "AmpliconDetector.findMatePrimer(rec: SamRecord)" should "find matches for the mate using the aligned/clipped coordinates" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=0, unclippedCoordinates=false)
    val builder  = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=100, cigar1="10S90M40H", start2=101, cigar2="100M")
    detector.findMatePrimer(r1).value shouldBe amplicon
    detector.findMatePrimer(r2).value shouldBe amplicon
  }

  it should "find matches for the mate using the unclipped coordinates" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=0, unclippedCoordinates=true)
    val builder  = new SamBuilder()
    val Seq(r1, r2) = builder.addPair(start1=110, cigar1="10S90M40H", start2=101, cigar2="90M10S")
    detector.findMatePrimer(r1).value shouldBe amplicon
    detector.findMatePrimer(r2).value shouldBe amplicon
  }

  "AmpliconDetector.find(r1: SamRecord, r2: SamRecord)" should "not matches for anything not mapped in FR pairs" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=5, unclippedCoordinates=false)
    val builder  = new SamBuilder()

    // not paired
    val frag = builder.addFrag(start=100, cigar="100M").value
    detector.find(r1=frag, r2=frag).isEmpty shouldBe true

    // one or both ends not mapped
    {
      val Seq(r1, r2) = builder.addPair(start1=100, cigar1="100M", start2=101, cigar2="100M")
      detector.find(r1=r1, r2=r2).value shouldBe amplicon

      // both unmapped
      r1.unmapped = true
      r2.unmapped = true
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true

      // r2 unmapped
      r1.mapped   = true
      r2.unmapped = true
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true

      // r1 unmapped
      r1.unmapped = true
      r2.mapped   = true
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }

    // different contigs
    {
      val Seq(r1, r2) = builder.addPair(start1=100, cigar1="100M", start2=101, cigar2="100M")
      detector.find(r1=r1, r2=r2).value shouldBe amplicon
      r1.refName = "chr1"
      r2.refName = "chr2"
      SamPairUtil.setMateInfo(r1.asSam, r2.asSam, true)
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }

    // not FR
    {
      val Seq(r1, r2) = builder.addPair(start1=100, cigar1="100M", start2=101, cigar2="100M", strand1=SamBuilder.Plus, strand2=SamBuilder.Plus)
      r1.isFrPair shouldBe false
      r2.isFrPair shouldBe false
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }
  }

  it should "find matches with aligned/clipped coordinates" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=5, unclippedCoordinates=false)
    val builder  = new SamBuilder()

    // exactly
    {
      val Seq(r1, r2) = builder.addPair(start1=100, cigar1="10H100M10H", start2=101, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicon
    }

    // slop
    {
      val Seq(r1, r2) = builder.addPair(start1=95, cigar1="10H100M10H", start2=106, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicon
    }

    // outside slop
    {
      val Seq(r1, r2) = builder.addPair(start1=94, cigar1="10H100M10H", start2=101, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }
    {
      val Seq(r1, r2) = builder.addPair(start1=100, cigar1="10H100M10H", start2=107, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }
  }

  it should "find matches with unclipped coordinates" in {
    val amplicon = Amplicon("chr1", 100, 120, 180, 200)
    val detector = AmpliconDetector(amplicons=Seq(amplicon), slop=5, unclippedCoordinates=true)
    val builder  = new SamBuilder()

    // exactly
    {
      val Seq(r1, r2) = builder.addPair(start1=110, cigar1="10H100M10H", start2=91, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicon
    }

    // slop
    {
      val Seq(r1, r2) = builder.addPair(start1=105, cigar1="10H100M10H", start2=96, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicon
    }

    // outside slop
    {
      val Seq(r1, r2) = builder.addPair(start1=104, cigar1="10H100M10H", start2=91, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }
    {
      val Seq(r1, r2) = builder.addPair(start1=110, cigar1="10H100M10H", start2=97, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).isEmpty shouldBe true
    }
  }

  it should "only return the hit with the smallest combined distance for multiple overlaps" in {
    val amplicons = IndexedSeq(Amplicon("chr1", 100, 120, 180, 200), Amplicon("chr1", 110, 130, 190, 210))
    val detector  = AmpliconDetector(amplicons=amplicons, slop=100, unclippedCoordinates=false)
    val builder   = new SamBuilder()

    // exactly to amplicon #1
    {
      val Seq(r1, r2) = builder.addPair(start1=100, cigar1="10H100M10H", start2=101, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicons.head
    }

    // one closer combined to amplicon #1
    {
      val Seq(r1, r2) = builder.addPair(start1=104, cigar1="10H100M10H", start2=105, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicons.head
    }

    // equal combined distance
    {
      val Seq(r1, r2) = builder.addPair(start1=105, cigar1="10H100M10H", start2=106, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).nonEmpty shouldBe true
    }

    // one closer to combined to amplicon #2
    {
      val Seq(r1, r2) = builder.addPair(start1=106, cigar1="10H100M10H", start2=107, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicons.last
    }

    // unbalanced
    {
      // read left primer is 4bp away from amplicon #1 left primer, read right primer is 17 bp from amplicon #1 right primer
      // read left primer is 6bp away from amplicon #2 left primer, read right primer is 7 bp from amplicon #2 right primer
      // so (4 + 17) = 21 > 13 = (6 + 7), so map to amplicon #2
      val Seq(r1, r2) = builder.addPair(start1=104, cigar1="10H100M10H", start2=118, cigar2="10H100M10H")
      detector.find(r1=r1, r2=r2).value shouldBe amplicons.last
    }
  }
}
