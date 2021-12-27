/*
 * The MIT License
 *
 * Copyright (c) 2019 Fulcrum Genomics LLC
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

package com.fulcrumgenomics.util

class UpdateGffContigNamesTest extends UpdateContigNamesSpec {

  private val gffLines =
    Seq(
      "##gff-version 3",
      "#!gff-spec-version 1.21",
      "#!processor NCBI annotwriter",
      "#!genome-build GRCh37.p13",
      "#!genome-build-accession NCBI_Assembly:GCF_000001405.25",
      "#!annotation-date",
      "#!annotation-source",
      "##sequence-region NC_000001.10 1 249250621",
      "##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=9606",
      "NC_000001.10\tRefSeq\tregion\t1\t249250621\t.\t+\t.\tID=id0;Dbxref=taxon:9606;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA",
      "NC_000001.10\tBestRefSeq\tgene\t11874\t14409\t.\t+\t.\tID=gene0;Dbxref=GeneID:100287102,HGNC:HGNC:37102;Name=DDX11L1;description=DEAD/H-box helicase 11 like 1;gbkey=Gene;gene=DDX11L1;gene_biotype=misc_RNA;pseudo=true",
      "NC_000001.10\tBestRefSeq\ttranscript\t11874\t14409\t.\t+\t.\tID=rna0;Parent=gene0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;Name=NR_046018.2;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2",
      "NC_000001.10\tBestRefSeq\texon\t11874\t12227\t.\t+\t.\tID=id1;Parent=rna0;Dbxref=GeneID:100287102,Genbank:NR_046018.2,HGNC:HGNC:37102;gbkey=misc_RNA;gene=DDX11L1;product=DEAD/H-box helicase 11 like 1;transcript_id=NR_046018.2",
      "##sequence-region NC_000002.10 1 243199373",
      "##sequence-region NC_000003.10 1 198022430",
      "##sequence-region NC_000004.10 1 191154276",
      "NC_000004.10\tRefSeq\tregion\t1\t191154276\t.\t+\t.\tID=id135365;Dbxref=taxon:9606;Name=4;chromosome=4;gbkey=Src;genome=chromosome;mol_type=genomic DNA",
      "NC_000004.10\tBestRefSeq\tgene\t53179\t88099\t.\t+\t.\tID=gene9164;Dbxref=GeneID:152687,HGNC:HGNC:27196,HPRD:15868;Name=ZNF595;description=zinc finger protein 595;gbkey=Gene;gene=ZNF595;gene_biotype=protein_coding",
      "NC_000004.10\tBestRefSeq\tmRNA\t53179\t88099\t.\t+\t.\tID=rna13646;Parent=gene9164;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Name=NM_001286054.1;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1",
      "NC_000004.10\tBestRefSeq\texon\t53179\t53385\t.\t+\t.\tID=id135366;Parent=rna13646;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1",
      "NC_000004.10\tBestRefSeq\texon\t59323\t59449\t.\t+\t.\tID=id135367;Parent=rna13646;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1",
      "NC_000004.10\tBestRefSeq\texon\t59951\t60046\t.\t+\t.\tID=id135368;Parent=rna13646;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1",
      "NC_000004.10\tBestRefSeq\texon\t83180\t83280\t.\t+\t.\tID=id135369;Parent=rna13646;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1",
      "NC_000004.10\tBestRefSeq\texon\t85622\t85995\t.\t+\t.\tID=id135370;Parent=rna13646;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1",
      "NC_000004.10\tBestRefSeq\texon\t85998\t86015\t.\t+\t.\tID=id135371;Parent=rna13646;Dbxref=GeneID:152687,Genbank:NM_001286054.1,HGNC:HGNC:27196,HPRD:15868;Note=The RefSeq transcript has 11 substitutions%2C 5 frameshifts%2C 1 non-frameshifting indel compared to this genomic sequence;exception=annotated by transcript or proteomic data;gbkey=mRNA;gene=ZNF595;product=zinc finger protein 595%2C transcript variant 4;transcript_id=NM_001286054.1"
    )

  private def gff = {
    val out = makeTempFile("test.", ".gff")
    Io.writeLines(out, gffLines)
    out
  }

  "UpdateGffContigNames" should "update a GFF" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateGffContigNames(
      input   = gff,
      dict    = pathToSequenceDictionary(),
      output  = output
    )

    executeFgbioTool(tool)

    val lines = Io.readLines(output).toSeq
    lines.filter(_.startsWith("##sequence-region")) should contain theSameElementsInOrderAs Seq(
      "##sequence-region chr1 1 249250621",
      "##sequence-region chr2 1 243199373",
      "##sequence-region chr3 1 198022430",
      "##sequence-region chr4 1 191154276"
    )
    lines.filterNot(_.startsWith("#")).flatMap(_.split('\t').headOption).distinct should contain theSameElementsInOrderAs Seq("chr1", "chr4")
  }

  it should "throw an exception if there are missing source contigs" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateGffContigNames(
      input   = gff,
      dict    = pathToSequenceDictionary(skipLast = true),
      output  = output
    )

    val ex = intercept[Exception] {executeFgbioTool(tool) }
    ex.getMessage should include ("Did not find contig")
  }

  it should "skip missing source contigs when using --skip-missing" in {
    val output = makeTempFile("test.", ".gff")
    val tool = new UpdateGffContigNames(
      input       = gff,
      dict        = pathToSequenceDictionary(skipLast = true),
      output      = output,
      skipMissing = true
    )

    executeFgbioTool(tool)

    val lines = Io.readLines(output).toSeq
    lines.filter(_.startsWith("##sequence-region")) should contain theSameElementsInOrderAs Seq(
      "##sequence-region chr1 1 249250621",
      "##sequence-region chr2 1 243199373",
      "##sequence-region chr3 1 198022430"
    )
    lines.filterNot(_.startsWith("#")).flatMap(_.split('\t').headOption).distinct should contain theSameElementsInOrderAs Seq("chr1")
  }
}