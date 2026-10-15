/**
 * @file queue
 * @author Paul Pillot <paul.pillot@libmol.org>
 * @license MIT
 * @copyright 2021
 * @description Queue JS implementation
 *
 */

/**
 * Double Ended Queue
 *
 * @export
 * @class DEQueue
 * @example let queue = new DEQueue(10);
 * queue.pushHead(1);
 * queue.pushHead(2);
 * let lTail = queue.popTail(); // lTail == 1
 * queue.size   // 1
 */
export class DEQueue<T> {
  store: T[];
  head: number;
  tail: number;
  storeSize: number;

  constructor(size: number) {
    this.storeSize = size;
    this.store = new Array(size);
    this.head = 0;
    this.tail = 0;
  }

  get size(): number {
    return (this.tail - this.head + this.storeSize) % this.storeSize;
  }

  get isEmpty(): boolean {
    return this.head === this.tail;
  }

  public popHead(): T | null {
    if (this.isEmpty) return null;

    const val = this.store[this.head];

    this.head++;
    if (this.head >= this.storeSize) this.head = 0;

    return val;
  }

  public popTail(): T | null {
    if (this.isEmpty) return null;

    this.tail--;
    if (this.tail < 0) this.tail = this.storeSize - 1;

    const val = this.store[this.tail];

    return val;
  }

  public pushHead(val: T) {
    this.head--;
    if (this.head < 0) this.head = this.storeSize - 1;

    this.store[this.head] = val;

    // handle overflow by shifting tail
    if (this.head === this.tail) {
      this.tail = (this.tail - 1 + this.storeSize) % this.storeSize;
    }
  }

  public pushTail(val: T) {
    this.store[this.tail] = val;

    this.tail++;
    if (this.tail >= this.storeSize) this.tail = 0;

    //handle overflow by shifting head
    if (this.head === this.tail) {
      this.head = (this.head - 1 + this.storeSize) % this.storeSize;
    }
  }

  public getHead(pos = 0): T {
    const lIdx = pos === 0 ? this.head : (this.head + pos + this.storeSize) % this.storeSize;
    return this.store[lIdx];
  }

  public getTail(pos = 0): T {
    const lIdx = (this.tail - 1 - pos + this.storeSize) % this.storeSize;
    return this.store[lIdx];
  }
}
