/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Scratch buffers which grew beyond this many entries (1 MB of 32-bit integers) are dropped after use, so that one huge call does not pin memory forever. */
export const MaxRetainedBufferCapacity = 1 << 18;

const Releasers = new Map<string, Set<() => void>>();

/** The scratch memory which code keeps between calls to avoid allocating (module level buffers): it is small after each call, and
 * `ScratchMemory.release(...contexts)` frees what the contexts kept (e.g. after loading a large file, when it is not needed for a while).
 * A context is the name of what uses the memory, e.g. `'schema-selection'`, the resolution of the rows of `StructureElement.Schema`. */
export const ScratchMemory = {
  /** Register a function which frees some scratch memory of a context (called by `release`). Return a function which unregisters it. */
  register(context: string, release: () => void): () => void {
    let releasers = Releasers.get(context);
    if (!releasers) Releasers.set(context, (releasers = new Set()));
    releasers.add(release);
    return () => {
      Releasers.get(context)?.delete(release);
    };
  },
  /** Free the scratch memory of the contexts, or of every context if there are none. It is allocated again when it is needed. */
  release(...contexts: string[]) {
    if (contexts.length > 0) {
      for (const context of contexts) for (const release of Releasers.get(context) ?? []) release();
    } else {
      for (const releasers of Releasers.values()) for (const release of releasers) release();
    }
  },
};

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
