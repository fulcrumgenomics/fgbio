/*
 * The MIT License
 *
 * Copyright (c) 2016 Fulcrum Genomics
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

import com.fulcrumgenomics.testing.UnitSpec

import scala.collection.mutable.ListBuffer

/**
  * Tests for CloseableBlockingQueue.
  */
class CloseableBlockingQueueTest extends UnitSpec {
  import QueueOp._

  object QueueOp {
    /** Trait that has a method that assess if the given queue has blocked threads. */
    trait RunnableWithHasQueuedThreads extends Runnable {
      protected def hasQueuedThreads[E](queue: CloseableBlockingQueue[E]): Boolean = {
        queue.lock.lock()
        val has = queue.lock.hasQueuedThreads
        queue.lock.unlock()
        has
      }
    }
    /** Runnable that waits until the queue has blocked threads, checks the given condition, then drains the queue. */
    class DrainQueue[E](queue: CloseableBlockingQueue[E], cond: () => Unit) extends RunnableWithHasQueuedThreads {
      def run(): Unit = {
        // wait until another thread is blocked on the queue
        while (!hasQueuedThreads(queue=queue)) { Thread.sleep(1) }
        // check the condition
        cond()
        // drain the queue
        while (queue.nonEmpty) queue.takeToMax(new ListBuffer[E]())
      }
    }
    /** Runnable that waits until the queue has blocked threads, checks the given condition, then fills the queue with the given element. */
    class FillQueue[E](queue: CloseableBlockingQueue[E], e: E, cond: () => Unit) extends RunnableWithHasQueuedThreads {
      def run(): Unit = {
        while (!hasQueuedThreads(queue=queue)) { Thread.sleep(1) }
        // check the condition
        cond()
        // fill the queue
        while (queue.remainingCapacity > 0) queue.putAll(Seq.range(0, queue.remainingCapacity, 1).map(i => e).toList)
      }
    }
    /** Runnable that applies the given function when run. */
    class DoQueueOp[E](queue: CloseableBlockingQueue[E], queueOp: CloseableBlockingQueue[E] => _) extends Runnable {
      def run(): Unit = queueOp(queue)
    }
    /** Tries to join the thread, and it fails, interrups. */
    private def joinOrInterrupt(thread: Thread, timeout: Int = 10): Unit = {
      thread.join(timeout)
      if (thread.isAlive) {
        thread.interrupt()
        thread.join()
      }
    }
    /** Runs threads where the blocker thread will initially block, and the unblocker will unblock the blocker thread. */
    def runBlockerAndUnblocker(blocker: Thread, unblocker: Thread): Unit = {
      blocker.start()
      unblocker.start()
      unblocker.join()
      joinOrInterrupt(blocker)
    }
  }

  "CloseableBlockingQueue" should "return the remaining capacity, size, and if it is empty or not" in {
    val queue = new CloseableBlockingQueue[Int](10)

    queue.remainingCapacity shouldBe 10
    queue.size() shouldBe 0
    queue.isEmpty shouldBe true
    queue.nonEmpty shouldBe false

    queue.put(1)
    queue.remainingCapacity shouldBe 9
    queue.size() shouldBe 1

    queue.putAll(Seq(1,2,3))
    queue.remainingCapacity shouldBe 6
    queue.size() shouldBe 4

    queue.take()
    queue.remainingCapacity shouldBe 7
    queue.size() shouldBe 3
    queue.isEmpty shouldBe false
    queue.nonEmpty shouldBe true
  }

  it should "put elements in the queue" in {
    val queue = new CloseableBlockingQueue[Int](10)
    val buffer = new ListBuffer[Int]()

    queue.put(1) shouldBe true
    queue.takeToMax(buffer) shouldBe 1
    buffer should contain theSameElementsInOrderAs List(1)
    buffer.clear()

    queue.putAll(Seq(1,2,3)) shouldBe true
    queue.takeToMax(buffer) shouldBe 3
    buffer should contain theSameElementsInOrderAs List(1,2,3)
    buffer.clear()

    queue.putSome(Seq(1,2,3)) shouldBe 3
    queue.takeToMax(buffer) shouldBe 3
    buffer should contain theSameElementsInOrderAs List(1,2,3)
    buffer.clear()
  }

  it should "take elements from the queue" in {
    val queue = new CloseableBlockingQueue[Int](10)
    val buffer = new ListBuffer[Int]()

    queue.put(1) shouldBe true
    queue.take() shouldBe 'defined

    queue.putAll(Seq(1,2,3)) shouldBe true
    queue.takeToMax(buffer) shouldBe 3
    buffer should contain theSameElementsInOrderAs List(1,2,3)
    buffer.clear()

    queue.putAll(Seq(1,2,3)) shouldBe true
    queue.takeToMax(buffer, max=2) shouldBe 2
    buffer should contain theSameElementsInOrderAs List(1,2)
    buffer.clear()
    queue.take() shouldBe 'defined
    queue.isEmpty shouldBe true

    queue.putAll(Seq(1,2,3)) shouldBe true
    queue.takeAtLeast(buffer, min=2) shouldBe 3
    buffer should contain theSameElementsInOrderAs List(1,2, 3)
    buffer.clear()
    queue.isEmpty shouldBe true
  }

  it should "not put or take elements after it is closed" in {
    val queue = new CloseableBlockingQueue[Int](10)
    val buffer = new ListBuffer[Int]()

    queue.putAll(Seq(1,2,3))
    queue.close()

    queue.put(1) shouldBe false
    queue.size() shouldBe 3

    queue.putAll(Seq(1,2,3)) shouldBe false
    queue.size() shouldBe 3

    queue.putSome(Seq(1,2,3)) shouldBe 0
    queue.size() shouldBe 3

    queue.take() shouldBe 'empty
    queue.size() shouldBe 3

    queue.takeToMax(buffer) shouldBe 0
    queue.size() shouldBe 3

    queue.takeAtLeast(buffer) shouldBe 0
    queue.size() shouldBe 3
  }

  "CloseableBlockingQueue.put()" should "block if there is no more capacity" in {
    val queue = new CloseableBlockingQueue[Int](1)
    queue.put(1)

    var putSuccess = false
    val queueOp = (queue: CloseableBlockingQueue[Int]) => { putSuccess = queue.put(2) }
    val put     = new Thread(new DoQueueOp[Int](queue=queue, queueOp=queueOp))
    val drain   = new Thread(new DrainQueue[Int](queue=queue, cond=() => putSuccess shouldBe false))

    runBlockerAndUnblocker(blocker=put, unblocker=drain)
    putSuccess shouldBe true
    queue.toSeq should contain theSameElementsInOrderAs Seq(2)
  }

  "CloseableBlockingQueue.putAll()" should "block when putting more elements into the queue than the queue's remaining capacity" in {
    val queue = new CloseableBlockingQueue[Int](2)
    queue.put(1)

    var putAllSuccess = false
    val queueOp = (queue: CloseableBlockingQueue[Int]) => { putAllSuccess = queue.putAll(Seq(2,3)) }
    val putAll  = new Thread(new DoQueueOp[Int](queue=queue, queueOp=queueOp))
    val drain   = new Thread(new DrainQueue[Int](queue=queue, cond=() => putAllSuccess shouldBe false))

    runBlockerAndUnblocker(blocker=putAll, unblocker=drain)
    putAllSuccess shouldBe true
    queue.toSeq should contain theSameElementsInOrderAs Seq(3)
  }

  "CloseableBlockingQueue.putSome()" should "block when the minimum to put is greater than the queue's remaining capacity" in {
    val queue = new CloseableBlockingQueue[Int](2)
    queue.put(1)

    var putSomeN = 0
    val queueOp = (queue: CloseableBlockingQueue[Int]) => { putSomeN = queue.putSome(Seq(2,3), min=2) }
    val putSome = new Thread(new DoQueueOp[Int](queue=queue, queueOp=queueOp))
    val drain   = new Thread(new DrainQueue[Int](queue=queue, cond=() => putSomeN shouldBe 0))

    runBlockerAndUnblocker(blocker=putSome, unblocker=drain)
    putSomeN shouldBe 2
    queue.size() shouldBe 1
    queue.toSeq should contain theSameElementsInOrderAs Seq(3)
  }

  it should "put fewer than requested but at least the minimum when the request to put is greater than the queue's remaining capacity" in {
    val queue = new CloseableBlockingQueue[Int](2)
    queue.put(1)
    queue.putSome(Seq(2, 3)) shouldBe 1
    queue.toSeq should contain theSameElementsInOrderAs Seq(1, 2)
  }

  "CloseableBlockingQueue.take()" should "block when there are no elements in the queue" in {
    val queue = new CloseableBlockingQueue[Int](1)

    var takeElement: Option[Int] = None
    val queueOp = (queue: CloseableBlockingQueue[Int]) => { takeElement = queue.take() }
    val take    = new Thread(new DoQueueOp[Int](queue=queue, queueOp=queueOp))
    val fill    = new Thread(new FillQueue[Int](queue=queue, e=1, cond=() => takeElement shouldBe 'empty))

    runBlockerAndUnblocker(blocker=take, unblocker=fill)
    takeElement shouldBe 'defined
    takeElement.get shouldBe 1
  }

  "CloseableBlockingQueue.takeToMax()" should "block when there are no elements in the queue" in {
    val queue = new CloseableBlockingQueue[Int](1)

    val buffer    = new ListBuffer[Int]()
    var numTaken  = 0
    val queueOp   = (queue: CloseableBlockingQueue[Int]) => { numTaken = queue.takeToMax(buffer) }
    val takeToMax = new Thread(new DoQueueOp[Int](queue=queue, queueOp=queueOp))
    val fill      = new Thread(new FillQueue[Int](queue=queue, e=1, cond=() => { numTaken shouldBe 0; buffer.isEmpty shouldBe true }))

    runBlockerAndUnblocker(blocker=takeToMax, unblocker=fill)
    numTaken shouldBe 1
    buffer should contain theSameElementsInOrderAs Seq(1)
  }

  it should "return only the maximum number, even if more elements are available" in {
    val queue = new CloseableBlockingQueue[Int](2)
    queue.putAll(Seq(1, 2))

    val buffer = new ListBuffer[Int]()
    queue.takeToMax(buffer, max=1) shouldBe 1
    buffer should contain theSameElementsInOrderAs Seq(1)
  }

  "CloseableBlockingQueue.takeAtLeast()" should "block when there are no elements in the queue" in {
    val queue = new CloseableBlockingQueue[Int](2)
    queue.put(0)

    val buffer      = new ListBuffer[Int]()
    var numTaken    = 0
    val queueOp     = (queue: CloseableBlockingQueue[Int]) => { numTaken = queue.takeAtLeast(buffer, min=2) }
    val takeAtLeast = new Thread(new DoQueueOp[Int](queue=queue, queueOp=queueOp))
    val fill        = new Thread(new FillQueue[Int](queue=queue, e=1, cond=() => { numTaken shouldBe 0; buffer.isEmpty shouldBe false; buffer.length shouldBe 1 }))

    runBlockerAndUnblocker(blocker=takeAtLeast, unblocker=fill)
    numTaken shouldBe 3
    buffer should contain theSameElementsInOrderAs Seq(0, 1, 1)
  }
}

