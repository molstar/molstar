/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Items grouped by (color, alpha) with a counting sort over the batch's distinct keys
// (seqspec:viewer#the-render-backend): the order the Canvas2D backend draws a batch in, shared by the backends that
// write a group at a time (Canvas2D's and the SVG's). O(items + groups), no allocation once the scratch has grown.

import { NO_COLOR } from "./backend.js";

export class ColorGroups {
  /** The items group by group: group q is `order[start[q] .. start[q + 1])`. */
  order = new Int32Array(256);
  start = new Int32Array(257);
  // The scratch: per item its group; per group its key; a stamped open-addressing table of keys.
  gid = new Int32Array(256);
  keyC = new Int32Array(256);
  keyK = new Int32Array(256);
  slots = new Int32Array(512);
  stamp = new Int32Array(512);
  gen = 0;

  /** The group count; NO_COLOR items are left out. Groups are in order of first appearance. */
  group(color: Int32Array, alpha: Uint8Array, n: number): number {
    if (this.gid.length < n) {
      // A power of two: the hash table masks its probes with `length - 1`.
      const m = 1 << (32 - Math.clz32(Math.max(n, this.gid.length * 2) - 1));
      this.gid = new Int32Array(m);
      this.keyC = new Int32Array(m);
      this.keyK = new Int32Array(m);
      this.start = new Int32Array(m + 1);
      this.order = new Int32Array(m);
      this.slots = new Int32Array(m * 2);
      this.stamp = new Int32Array(m * 2);
    }
    const { gid, keyC, keyK, start, slots, stamp, order } = this;
    const mask = slots.length - 1;
    if (++this.gen === 0x7fffffff) {
      stamp.fill(0);
      this.gen = 1;
    }
    const gen = this.gen;
    let groups = 0;
    let lc = 0;
    let lk = -1;
    let lg = -1;
    for (let i = 0; i < n; i++) {
      const c = color[i];
      if (c === NO_COLOR) {
        gid[i] = -1;
        continue;
      }
      const k = alpha[i];
      if (c !== lc || k !== lk) {
        let h = (Math.imul(c, 0x9e3779b1) ^ k) & mask;
        while (stamp[h] === gen && (keyC[slots[h]] !== c || keyK[slots[h]] !== k)) h = (h + 1) & mask;
        if (stamp[h] !== gen) {
          stamp[h] = gen;
          slots[h] = groups;
          keyC[groups] = c;
          keyK[groups] = k;
          start[groups++] = 0;
        }
        lg = slots[h];
        lc = c;
        lk = k;
      }
      gid[i] = lg;
      start[lg]++;
    }
    // Counts to starts, then each item to its group's next place; the starts are restored by the shift.
    let sum = 0;
    for (let q = 0; q < groups; q++) {
      const m = start[q];
      start[q] = sum;
      sum += m;
    }
    start[groups] = sum;
    for (let i = 0; i < n; i++) if (gid[i] >= 0) order[start[gid[i]]++] = i;
    for (let q = groups; q > 0; q--) start[q] = start[q - 1];
    start[0] = 0;
    return groups;
  }
}
