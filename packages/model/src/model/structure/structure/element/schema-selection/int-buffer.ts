/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Scratch buffers which grew beyond this many entries are dropped after use, so that one huge call does not pin memory forever. */
export const MaxRetainedBufferCapacity = 1 << 22;

/** Growable buffer of 32-bit integers with an explicit element count, for reuse across calls.
 *
 * (Reusing `array.length = 0` does not work for this: V8 releases the backing store when an array is truncated, so
 * every refill allocates again. A buffer which only resets its counter keeps its capacity.) */
export class IntBuffer {
  data: Int32Array;
  n = 0;

  constructor(capacity = 256) {
    this.data = new Int32Array(capacity);
  }

  clear() {
    this.n = 0;
  }

  /** Make sure `capacity` values fit without growing (keeps the current contents). */
  ensureCapacity(capacity: number) {
    if (capacity <= this.data.length) return;
    const data = new Int32Array(Math.max(capacity, this.data.length * 2));
    data.set(this.data.subarray(0, this.n));
    this.data = data;
  }

  push(value: number) {
    if (this.n === this.data.length) this.ensureCapacity(this.n + 1);
    this.data[this.n++] = value;
  }

  /** The filled part as a typed array view (shares memory with the buffer, valid until the buffer is modified). */
  view(): Int32Array {
    return this.data.subarray(0, this.n);
  }

  /** Free the memory if the capacity grew above `limit` (so one huge call does not pin memory forever). */
  trim(limit: number) {
    if (this.data.length > limit) {
      this.data = new Int32Array(256);
      this.n = 0;
    }
  }
}
