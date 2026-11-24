/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Hits held columnar (_spec/README.md "Search", Hits): `Int32Array` columns grown by doubling, `Int8Array` strands; no
// object a hit.

import type { SearchHits } from './types.js';

/** A growing list of hits: what a batch, and the viewer's result, hold. */
export class HitList implements SearchHits {
  length = 0;
  target: Int32Array;
  start: Int32Array;
  end: Int32Array;
  strand: Int8Array;

  constructor(capacity = 16) {
    const c = Math.max(1, capacity);
    this.target = new Int32Array(c);
    this.start = new Int32Array(c);
    this.end = new Int32Array(c);
    this.strand = new Int8Array(c);
  }

  push(target: number, start: number, end: number, strand: number): void {
    if (this.length === this.start.length) this.grow(this.length * 2);
    const i = this.length++;
    this.target[i] = target;
    this.start[i] = start;
    this.end[i] = end;
    this.strand[i] = strand;
  }

  /** Appends `count` hits of another list from `from`. */
  append(other: SearchHits, from = 0, count = other.length - from): void {
    if (count <= 0) return;
    if (this.length + count > this.start.length) this.grow(Math.max(this.length + count, this.length * 2));
    const at = this.length;
    this.target.set(other.target.subarray(from, from + count), at);
    this.start.set(other.start.subarray(from, from + count), at);
    this.end.set(other.end.subarray(from, from + count), at);
    this.strand.set(other.strand.subarray(from, from + count), at);
    this.length += count;
  }

  clear(): void {
    this.length = 0;
  }

  /** A list of exactly these hits (columns trimmed to the length). */
  trimmed(): SearchHits {
    const n = this.length;
    return {
      length: n,
      target: this.target.slice(0, n),
      start: this.start.slice(0, n),
      end: this.end.slice(0, n),
      strand: this.strand.slice(0, n),
    };
  }

  private grow(capacity: number): void {
    const t = new Int32Array(capacity);
    const s = new Int32Array(capacity);
    const e = new Int32Array(capacity);
    const d = new Int8Array(capacity);
    t.set(this.target.subarray(0, this.length));
    s.set(this.start.subarray(0, this.length));
    e.set(this.end.subarray(0, this.length));
    d.set(this.strand.subarray(0, this.length));
    this.target = t;
    this.start = s;
    this.end = e;
    this.strand = d;
  }
}

/** No hits. */
export const NoHits: SearchHits = {
  length: 0,
  target: new Int32Array(0),
  start: new Int32Array(0),
  end: new Int32Array(0),
  strand: new Int8Array(0),
};
