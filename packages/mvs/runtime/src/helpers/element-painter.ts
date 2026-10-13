/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** Shared "next free element" buffer, only used during a painting pass and never returned to callers. */
let Next = new Int32Array(0);
/** Buffers larger than this are dropped after use so that one huge call does not pin memory forever. */
const MaxRetainedCapacity = 1 << 22;

/** Assigns a value to each element of `[0, n)`, where every element gets the value of the *first* painted range which contains it
 * (to give the *last* row of an annotation priority, paint the rows in reverse order).
 *
 * Painted elements are skipped through a union-find array of "next unpainted element" (with path compression), so every
 * element is written at most once however many ranges overlap it: painting is O(n + number of ranges) instead of
 * O(sum of range sizes). `isComplete` tells when all elements are painted and the remaining ranges cannot change anything. */
export class ElementPainter {
  private out: Int32Array | null = null;
  private n = 0;
  private remaining = 0;

  /** Start painting `n` elements (forgets the previous result). */
  begin(n: number) {
    this.n = n;
    this.remaining = n;
    this.out = null;
  }

  /** True when every element has a value (or there are no elements), so painting more ranges is pointless. */
  get isComplete() {
    return this.remaining === 0;
  }

  /** Paint the first `count` ranges `[from[i], to[i])` with `value` (only the elements which are not painted yet). */
  paint(from: ArrayLike<number>, to: ArrayLike<number>, count: number, value: number) {
    if (count === 0 || this.remaining === 0) return;
    const n = this.n;
    if (!this.out) {
      this.out = new Int32Array(n).fill(-1);
      if (Next.length < n + 1) Next = new Int32Array(Math.max(n + 1, 1024));
      Next.fill(0, 0, n + 1); // 0 = not painted (a painted element `i` stores a larger index: the next candidate)
    }
    const out = this.out;
    const next = Next;
    for (let r = 0; r < count; r++) {
      const stop = to[r];
      let i = findFree(next, from[r]);
      while (i < stop) {
        out[i] = value;
        next[i] = i + 1;
        this.remaining--;
        i = findFree(next, i + 1);
      }
      if (this.remaining === 0) return;
    }
  }

  /** The painted values (`-1` for elements no range contained), or `null` if nothing was painted. */
  result(): Int32Array | null {
    return this.out;
  }

  /** Free the shared buffer if it grew large. */
  static release() {
    if (Next.length > MaxRetainedCapacity) Next = new Int32Array(0);
  }
}

/** First element `>= i` which is not painted yet (`n` if none), compressing the path on the way. */
function findFree(next: Int32Array, i: number): number {
  let root = i;
  while (next[root] !== 0) root = next[root];
  while (next[i] !== 0) {
    const following = next[i];
    next[i] = root;
    i = following;
  }
  return root;
}
