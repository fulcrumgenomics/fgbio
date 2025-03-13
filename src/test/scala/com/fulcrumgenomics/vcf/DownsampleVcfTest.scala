package com.fulcrumgenomics.vcf

import com.fulcrumgenomics.testing.VcfBuilder
import com.fulcrumgenomics.testing.VcfBuilder.Gt
import com.fulcrumgenomics.util.Metric
import com.fulcrumgenomics.vcf.api.Allele.SimpleAllele
import com.fulcrumgenomics.vcf.api.{Allele, AlleleSet, Genotype, Variant}
import com.fulcrumgenomics.testing.UnitSpec
import com.fulcrumgenomics.vcf.DownsampleVcf.{Likelihoods, computePls, downsampleAndRegenotype}

import scala.util.Random

class DownsampleVcfTest extends UnitSpec {

  private val dict = VcfBuilder(Seq.empty).header.dict
  
  private def winnow(variants: Iterator[Variant], window: Int = 10): Iterator[Variant] = {
    DownsampleVcf.winnowVariants(variants, window, dict)
  }
  
  "DownsampleVcf.winnowVariants" should "return an empty iterator when given an empty iterator" in {
    winnow(Iterator.empty, window = 0) shouldBe Symbol("empty")
  }

  it should "return a single variant when given a single variant" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","G"))
    winnow(builder.iterator, window=0).toSeq should contain theSameElementsInOrderAs builder.iterator.toSeq
  }

  it should "throw an exception if the window is negative" in {
    val builder = VcfBuilder(Seq.empty)
    val ex = intercept[Exception] { winnow(builder.iterator, window= -2).toSeq }
    ex.getMessage should include ("negative")
  }

  it should "keep both variants if distance is one more than the boundary" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","G"))
    builder.add(pos=14, alleles=Seq("G","C"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs builder.iterator.toSeq
  }

  it should "keep both variants if distance is the boundary" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","G"))
    builder.add(pos=13, alleles=Seq("G","C"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs builder.iterator.toSeq
  }

  it should "keep the first variant if distance is within specified window" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","G"))
    builder.add(pos=9, alleles=Seq("G","C"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs Iterator(builder.iterator.next()).toSeq
  }

  it should "keep both variants if genomic loci is in range, but on different chromosomes" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(chrom="chr1", pos=3, alleles=Seq("A","G"))
    builder.add(chrom="chr2", pos=9, alleles=Seq("G","C"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs builder.iterator.toSeq
  }

  it should "keep the first variant if three variants are within the specified window" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","G"))
    builder.add(pos=9, alleles=Seq("G","C"))
    builder.add(pos=12, alleles=Seq("G","T"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs Iterator(builder.iterator.next()).toSeq
  }

  it should "remove the middle variant if variants a and b are within window, b and c are within window, and a and c are far enough apart" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","G"))
    builder.add(pos=11, alleles=Seq("G","C"))
    builder.add(pos=14, alleles=Seq("G","T"))
    val answer = builder.iterator.filterNot(_.pos==11)
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs answer.toSeq
  }

  it should "throw an exception if the chromosomes are out of order" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(chrom="chr2", pos=3, alleles=Seq("A","G"))
    builder.add(chrom="chr1", pos=9, alleles=Seq("G","C"))
    val ex = intercept[Exception] { winnow(builder.iterator.toSeq.reverseIterator, window=10).toSeq }
    ex.getMessage should include ("variants out of order")
  }

  it should "throw an exception if the positions within a chromosome are out of order" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=9, alleles=Seq("A","G"))
    builder.add(pos=3, alleles=Seq("G","C"))
    val ex = intercept[Exception] { winnow(builder.iterator.toSeq.reverseIterator, window=10).toSeq }
    ex.getMessage should include ("variants out of order")
  }

  it should "keep the first variant if there are duplicate variants" in {
    // VcfBuilder doesn't support duplicate variants
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=9, alleles=Seq("A","G"))
    val builder2 = VcfBuilder(Seq.empty)
    builder2.add(pos=9, alleles=Seq("A","G"))
    val builder3 = Iterator(builder.iterator.next()) ++ Iterator(builder2.iterator.next())
    winnow(builder3.iterator, window=10).toSeq should contain theSameElementsInOrderAs Iterator(builder.iterator.next()).toSeq
  }

  it should "keep the first variant if there are variants at the same contig/position" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=9, alleles=Seq("A","G"))
    val builder2 = VcfBuilder(Seq.empty)
    builder2.add(pos=9, alleles=Seq("A","T"))
    val builder3 = Iterator(builder.iterator.next()) ++ Iterator(builder2.iterator.next())
    winnow(builder3.iterator, window=10).toSeq should contain theSameElementsInOrderAs Iterator(builder.iterator.next()).toSeq
  }

  it should "keep both variants if there is an insertion and the following SNP is outside the window" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","AGGG"))
    builder.add(pos=14, alleles=Seq("C","T"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs builder.toSeq
  }

  it should "keep the first variant if there is an insertion and the following SNP is within the window" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("A","AGGG"))
    builder.add(pos=12, alleles=Seq("C","T"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs Iterator(builder.iterator.next()).toSeq
  }

  it should "keep both variants if there is a deletion and the following SNP is outside the window" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("AGGG","A"))
    builder.add(pos=16, alleles=Seq("C","G"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs builder.toSeq
  }

  it should "keep the first variant if there is a deletion and the following SNP is within the window" in {
    val builder = VcfBuilder(Seq.empty)
    builder.add(pos=3, alleles=Seq("AGGG","A"))
    builder.add(pos=15, alleles=Seq("C","G"))
    winnow(builder.iterator, window=10).toSeq should contain theSameElementsInOrderAs Iterator(builder.iterator.next()).toSeq
  }

  /*
    Testing DownsampleVcf.downsampleADs
  */

  val random = new Random(60)

  "DownsampleVcf.downsampleADs" should "return the same # of ADs as was given" in {
    val ad = IndexedSeq(1, 2)
    val expected = DownsampleVcf.downsampleADs(ad, proportion = 0.1, random)
    expected.length shouldBe ad.length
  }

  it should "downsample to roughly 1% of the previous allele number" in {
    val ad = IndexedSeq(10000, 25000)
    val expected = DownsampleVcf.downsampleADs(ad, proportion = 0.1, random)
    expected should contain theSameElementsInOrderAs Seq(1003, 2523)
  }

  it should "return 0 if the input allele depths are 0" in {
    val ad = IndexedSeq(0, 0)
    val expected = DownsampleVcf.downsampleADs(ad, proportion = 0.1, random)
    expected should contain theSameElementsInOrderAs Seq(0,0)
  }

  it should "return a small number or 0 if the input ADs are small" in {
    val ad = IndexedSeq(1, 10)
    val expected = DownsampleVcf.downsampleADs(ad, proportion = 0.1, random)
    expected should contain theSameElementsInOrderAs Seq(0,1)
  }

  it should "return the same values if the proportion is 1" in {
    val ad = IndexedSeq(100, 100)
    val expected = DownsampleVcf.downsampleADs(ad, proportion = 1, random)
    expected should contain theSameElementsInOrderAs Seq(100, 100)
  }

  it should "throw an exception if the proportion is greater than 1" in {
    val ad = IndexedSeq(100, 200)
    val ex = intercept[Exception] { DownsampleVcf.downsampleADs(ad, proportion = 1.1, random) }
    ex.getMessage should include ("proportion must be less than 1")
  }

  "DownsampleVcf.computePls" should "return new PLs that are not always 0,0,0" in {
    val ads = IndexedSeq[Int](0, 100)
    val expected = IndexedSeq(1996, 301, 0)
    val newlikelihoods = computePls(ads)
    newlikelihoods should contain theSameElementsInOrderAs expected
  }

  /*
    Testing DownsampleVcf.Likelihoods
   */

  "DownsampleVcf.Likelihoods" should "return ref if all allele depths are zero" in {
    val likelihood = Likelihoods(alleleDepthA=0, alleleDepthB=0)
    val expected = IndexedSeq[Int](0, 0, 0)
    likelihood.pls.length shouldBe expected.length
    likelihood.pls should contain theSameElementsInOrderAs expected
  }

  it should "return a likelihood of 0 for AA if there are only ref alleles observed" in {
    val likelihood = Likelihoods(alleleDepthA = 10, alleleDepthB = 0)
    val expected = IndexedSeq[Int](0, 30, 200)
    likelihood.pls should contain theSameElementsInOrderAs expected
  }

  it should "return a likelihood of 0 for BB if there are only alt alleles observed" in {
    val likelihood = Likelihoods(alleleDepthA = 0, alleleDepthB = 10)
    val expected = IndexedSeq[Int](200, 30, 0)
    likelihood.pls should contain theSameElementsInOrderAs expected
  }

  it should "return a likelihood of 0 for AB if there are an equal number of ref and alt alleles" in {
    val likelihood = Likelihoods(alleleDepthA = 5, alleleDepthB = 5)
    val expected = IndexedSeq[Int](70, 0, 70)
    likelihood.pls should contain theSameElementsInOrderAs expected
  }

  it should "return a likelihood of 0 for AA if the AD A >> AD B" in {
    val likelihood = Likelihoods(alleleDepthA = 15, alleleDepthB = 2)
    likelihood.pls(0) == 0
  }

  it should "return a likelihood of 0 for AB if AD.A and AD.B are similar but not equal" in {
    val likelihood = Likelihoods(alleleDepthA = 15, alleleDepthB = 17)
    likelihood.pls(1) == 0
  }

  it should "return a likelihood of 0 for BB if AD.B >> AD.A but neither are 0" in {
    val likelihood = Likelihoods(alleleDepthA = 3, alleleDepthB = 30)
    likelihood.pls(2) == 0
  }

  it should "return correct values when there are very few reads" in {
    Likelihoods(0, 0).pls should contain theSameElementsInOrderAs IndexedSeq(0, 0, 0)
    Likelihoods(1, 0).pls should contain theSameElementsInOrderAs IndexedSeq(0, 3, 20)
    Likelihoods(1, 1).pls should contain theSameElementsInOrderAs IndexedSeq(14, 0, 14)
    Likelihoods(0, 2).pls should contain theSameElementsInOrderAs IndexedSeq(40, 6, 0)
    Likelihoods(1, 2).pls should contain theSameElementsInOrderAs IndexedSeq(31, 0, 11)
  }


  /*
  testing DownsampleVcf.downsampleAndRegenotype on downsampleAndRegenotypes
   */
  private def makeGt(ref: String, alt: String, ads: IndexedSeq[Int], sample: String ="test"): Genotype = {
    Genotype(alleles=AlleleSet(ref=SimpleAllele(ref), alts=IndexedSeq(Allele(alt))),
      sample=sample,
      calls=IndexedSeq[Allele](Allele(ref), Allele(alt)),
      attrs=Map("AD" -> ads, "PL" -> Likelihoods(alleleDepthA = ads(0), alleleDepthB = ads(1))))
  }

  
  "DownsampleVcf.downsampleAndRegneotype(Genotype)" should "return no call if all allele depths are zero" in {
    val geno = makeGt(ref="A", alt="T", ads=IndexedSeq(0,0))
    val newGeno = downsampleAndRegenotype(gt=geno, proportion=0.01, random = new Random(42), epsilon = 0.01)
    val expected = IndexedSeq(Allele("."), Allele("."))
    newGeno.calls should contain theSameElementsInOrderAs expected
  }

  it should "return two ref alleles if the ref AD is much larger than the alt AD" in {
    val geno = makeGt(ref="A", alt="T", ads=IndexedSeq(100,0))
    val newGeno = downsampleAndRegenotype(gt=geno, proportion=0.1, random = new Random(42), epsilon = 0.01)
    val expected = IndexedSeq(Allele("A"), Allele("A"))
    newGeno.calls should contain theSameElementsInOrderAs expected
  }

  it should "return two alt alleles if the alt AD is greater than the ref AD" in {
    val geno = makeGt(ref="A", alt="T", ads=IndexedSeq(0,100))
    val newGeno = downsampleAndRegenotype(gt=geno, proportion=0.1, random = new Random(42), epsilon = 0.01)
    val expected = IndexedSeq(Allele("T"), Allele("T"))
    newGeno.calls should contain theSameElementsInOrderAs expected
  }

  it should "return two alt alleles if ref and alt AD > 0 but ref << alt" in {
    val geno = makeGt(ref="A", alt="T", ads=IndexedSeq(30,200))
    val newGeno = downsampleAndRegenotype(gt=geno, proportion=0.1, random = new Random(42), epsilon = 0.01)
    val expected = IndexedSeq(Allele("T"), Allele("T"))
    newGeno.calls should contain theSameElementsInOrderAs expected
  }

  it should "return het if ref and alt ADs are similar but ref < alt" in {
    val geno = makeGt(ref="A", alt="T", ads=IndexedSeq(190,200))
    val newGeno = downsampleAndRegenotype(gt=geno, proportion=0.1, random = new Random(42), epsilon = 0.01)
    val expected = IndexedSeq(Allele("A"), Allele("T"))
    newGeno.calls should contain theSameElementsInOrderAs expected
  }


  it should "return a ref and alt allele if the ref and alt ADs are the same" in {
    val geno = makeGt(ref="A", alt="T", ads=IndexedSeq(100,100))
    val newGeno = downsampleAndRegenotype(gt=geno, proportion=0.1, random = new Random(42), epsilon = 0.01)
    val expected = IndexedSeq(Allele("A"), Allele("T"))
    newGeno.calls should contain theSameElementsInOrderAs expected
  }

  /*
    testing DownsampleVcf.downsampleAndRegenotype on Variant
   */

  private def makeVariant(ref: String, alt: String, sample: String = "test", ads: IndexedSeq[Int]): Variant = {
    Variant(chrom="1",
            pos=10,
            alleles=AlleleSet(ref=Allele(ref), alts=Allele(alt)),
            genotypes=Map(sample -> makeGt(ref=ref, alt=alt, ads=ads, sample = sample))
    )
  }

  "DownsampleVcf.downsampleAndRegenotype(Variant)" should "return no call alleles if depths are 0" in {
    val variant = makeVariant(ref="A", alt="T", ads=IndexedSeq(0,0))
    val newVariant = DownsampleVcf.downsampleAndRegenotype(variant=variant, proportions = Map("test" -> 0.1), random = new Random(42))
    val expected = IndexedSeq(Allele("."), Allele("."))
    newVariant.genotypes("test").calls should contain theSameElementsInOrderAs expected
  }

  it should "return ref alleles if ref AD > 0 and alt is 0" in {
    val variant = makeVariant(ref="A", alt="T", ads=IndexedSeq(100,0))
    val newVariant = DownsampleVcf.downsampleAndRegenotype(variant=variant, proportions = Map("test" -> 0.1), random = new Random(42))
    val expected = IndexedSeq(Allele("A"), Allele("A"))
    newVariant.genotypes("test").calls should contain theSameElementsInOrderAs expected
  }

  it should "return ref if ref and alt > 0 but ref > alt" in {
    val variant = makeVariant(ref="A", alt="T", ads=IndexedSeq(200,20))
    val newVariant = DownsampleVcf.downsampleAndRegenotype(variant=variant, proportions = Map("test" -> 0.1), random = new Random(42))
    val expected = IndexedSeq(Allele("A"), Allele("A"))
    newVariant.genotypes("test").calls should contain theSameElementsInOrderAs expected
  }

  it should "return alts if alt AD > 0 and ref AD = 0" in {
    val variant = makeVariant(ref="A", alt="T", ads=IndexedSeq(0,100))
    val newVariant = DownsampleVcf.downsampleAndRegenotype(variant=variant, proportions = Map("test" -> 0.1), random = new Random(42))
    val expected = IndexedSeq(Allele("T"), Allele("T"))
    newVariant.genotypes("test").calls should contain theSameElementsInOrderAs expected
  }

  it should "return het if ref and alt ADs are the same" in {
    val variant = makeVariant(ref="A", alt="T", ads=IndexedSeq(500,500))
    val newVariant = DownsampleVcf.downsampleAndRegenotype(variant=variant, proportions = Map("test" -> 0.1), random = new Random(42))
    val expected = IndexedSeq(Allele("A"), Allele("T"))
    newVariant.genotypes("test").calls should contain theSameElementsInOrderAs expected
  }

  private val sample = "test1"
  private val builder = VcfBuilder(samples=Seq(sample))
  builder.add(chrom="chr1", pos=100, id="1", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="0/1", ad=Seq(1000,1000), pl=Seq(13936,0,13936)))) // should stay
  builder.add(chrom="chr1", pos=200, id="2", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="0/0", ad=Seq(1000,0), pl=Seq(0,3010,19956)))) // should be removed
  builder.add(chrom="chr1", pos=400, id="3", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="1/1", ad=Seq(0,1000), pl=Seq(19956,3010,0)))) // should stay
  builder.add(chrom="chr1", pos=600, id="4", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="1/1", ad=Seq(0,3), pl=Seq(60,9,0)))) // should be removed
  builder.add(chrom="chr1", pos=800, id="5", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="1/1", ad=Seq(20,1000), pl=Seq(19557,2671,0)))) // should stay
  builder.add(chrom="chr1", pos=1000, id="6", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="1/1", ad=Seq(1000,10), pl=Seq(0,2841,19757)))) // should stay
  builder.add(chrom="chr1", pos=1200, id="7", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="1/1", ad=Seq(1,2), pl=Seq(31,0,11)))) // should be removed
  builder.add(chrom="chr1", pos=1360, id="8", alleles=Seq("A", "C"), info=Map(),
    gts=Seq(Gt(sample=sample, gt="1/1", ad=Seq(800,900), pl=Seq(12843,0,10848)))) // should stay
  val inVcf = builder.toTempFile()
  // Make the metadata file
  val metadata = makeTempFile("metadata", ".txt")
  Metric.write(metadata, Seq(Sample(SAMPLE_NAME = sample, BASE_COUNT = 100)))

  "DownsampleVcf" should "write a new vcf with downsampled genotypes when provided a vcf" in {
    List(true, false).foreach(
      use_metdata => {
        // Construct the input VCF
        val outVcf = makeTempFile("out", ".vcf.gz")
        if(use_metdata) {
          new DownsampleVcf(input=inVcf,
                            output=outVcf,
                            metadata=Some(metadata),
                            downsampleToBases=Some(1),
                            windowSize=150).execute()
        } else {
          new DownsampleVcf(input=inVcf,
                            output=outVcf,
                            originalBases=Some(100),
                            downsampleToBases=Some(1),
                            windowSize=150).execute()
        }

        val vs = readVcfRecs(outVcf)
        vs should have length 5

        val ad0 = vs(0).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad0(0) < 110 shouldBe true;
        ad0(1) < 110 shouldBe true
        val pl0 = vs(0).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl0(1) shouldBe 0

        val ad1 = vs(1).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad1(0) shouldBe 0;
        ad1(1) < 110 shouldBe true
        val pl1 = vs(1).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl1(2) shouldBe 0

        val ad2 = vs(2).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad2(0) < 30 shouldBe true;
        ad2(1) < 110 shouldBe true
        val pl2 = vs(2).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl2(2) shouldBe 0

        val ad3 = vs(3).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad3(0) < 110 shouldBe true;
        ad3(1) shouldBe 0;
        val pl3 = vs(3).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl3(0) shouldBe 0

        val ad4 = vs(4).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad4(0) < 100 shouldBe true;
        ad4(1) < 100 shouldBe true;
        ad4(0) > 1 shouldBe true;
        ad4(1) > 2 shouldBe true;
        val pl4 = vs(4).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl4(1) shouldBe 0
      }
    )

  }

  "DownsampleVcf" should "write a new vcf with downsampled genotypes when provided a vcf, keeping nocalls" in {
    // Construct the input VCF
    List(true, false).foreach(
      use_metdata => {
        // Construct the input VCF
        val outVcf = makeTempFile("out", ".vcf.gz")
        if (use_metdata) {
          new DownsampleVcf(input=inVcf, 
                            output=outVcf,
                            metadata=Some(metadata),
                            downsampleToBases=Some(1),
                            writeNoCall=true,
                            windowSize=150).execute()
        } else {
          new DownsampleVcf(input=inVcf,
                            output=outVcf,
                            originalBases=Some(100),
                            downsampleToBases=Some(1),
                            writeNoCall=true,
                            windowSize=150).execute()
        }

        val vs = readVcfRecs(outVcf)
        vs should have length 7

        val ad0 = vs(0).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad0(0) < 110 shouldBe true;
        ad0(1) < 110 shouldBe true
        val pl0 = vs(0).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl0(1) shouldBe 0

        val ad1 = vs(1).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad1(0) shouldBe 0;
        ad1(1) < 110 shouldBe true
        val pl1 = vs(1).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl1(2) shouldBe 0

        val ad2 = vs(2).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad2(0) shouldBe 0;
        ad2(1) shouldBe 0
        val pl2 = vs(2).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl2(0) shouldBe 0;
        pl2(1) shouldBe 0;
        pl2(2) shouldBe 0

        val ad3 = vs(3).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad3(0) < 30 shouldBe true;
        ad3(1) < 110 shouldBe true
        val pl3 = vs(3).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl3(2) shouldBe 0

        val ad4 = vs(4).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad4(0) < 110 shouldBe true;
        ad4(1) shouldBe 0;
        val pl4 = vs(4).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl4(0) shouldBe 0

        val ad5 = vs(5).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad5(0) shouldBe 0;
        ad5(1) shouldBe 0
        val pl5 = vs(5).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl5(0) shouldBe 0;
        pl5(1) shouldBe 0;
        pl5(2) shouldBe 0

        val ad6 = vs(6).genotypes("test1")[IndexedSeq[Int]]("AD")
        ad6(0) < 100 shouldBe true;
        ad6(1) < 100 shouldBe true;
        ad6(0) > 1 shouldBe true;
        ad6(1) > 2 shouldBe true;
        val pl6 = vs(6).genotypes("test1")[IndexedSeq[Int]]("PL")
        pl6(1) shouldBe 0
      }
    )
  }
}

