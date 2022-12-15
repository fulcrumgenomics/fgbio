package com.fulcrumgenomics.personal.yfarjoun

import com.fulcrumgenomics.FgBioDef.DirPath
import com.fulcrumgenomics.bam.Bams.templateIterator
import com.fulcrumgenomics.bam.{RandomizeBam, Template}
import com.fulcrumgenomics.bam.api.{SamOrder, SamRecord, SamRecordCodec}
import com.fulcrumgenomics.commons.collection.SelfClosingIterator
import com.fulcrumgenomics.commons.util.Logger
import com.fulcrumgenomics.util.Sorter
import htsjdk.samtools.SAMFileHeader
import htsjdk.samtools.util.Murmur3

object RandomTemplate {

  /** Return an iterator over records sorted and grouped into [[Template]] objects. Although a queryname sort is
    * guaranteed, the sort order may not be consistent with other queryname sorting implementations, especially in other
    * tool kits. See [[templateIterator]] for a [[Template]] iterator which emits templates in a non-guaranteed sort
    * order.
    *
    * @see [[templateIterator]]
    *
    * @param iterator an iterator from which to consume records
    * @param header the header associated with the records
    * @param maxInMemory the maximum number of records to keep and sort in memory, if sorting is needed
    * @param tmpDir a temp directory to use for temporary sorting files if sorting is needed
    * @return an Iterator of queryname sorted Template objects
    */
  def templateSortedIterator(iterator: Iterator[SamRecord],
                             header: SAMFileHeader,
                             maxInMemory: Int,
                             logger: Logger): SelfClosingIterator[Template] = {

    val sorter = RandomizeBam.Randomize(iterator, header, 42, true,logger)
    val sortedIterator = sorter.iterator

    val iter = new Iterator[Template] {
      override def hasNext: Boolean = sortedIterator.hasNext
      override def next(): Template = {
        require(hasNext, "next() called on empty iterator")
        Template(sortedIterator)
      }
    }

    new SelfClosingIterator[Template](iter, ()=>{sorter.close()})
  }
}
