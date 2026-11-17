/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Lanes and the hit index of marks in lanes (_spec/README.md "Render", Lanes and Hits).
//
// Sticky lanes: an item keeps its lane while it is in view, so items never jump when a neighboring region arrives; a
// new item takes the first lane free over its extent; the lanes are repacked only when the view has settled and a
// fresh packing needs fewer lanes than the kept ones (they became sparse).

type Key = string | number;

/** Whether `[s, e)` overlaps a lane's intervals (`starts` sorted; intervals in a lane never overlap). */
function overlaps(starts: number[], ends: number[], s: number, e: number): boolean {
  let lo = 0;
  let hi = starts.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (starts[mid] < e) lo = mid + 1;
    else hi = mid;
  }
  // the interval starting last before `e` is the only one that can reach past `s`
  return lo > 0 && ends[lo - 1] > s;
}

function insert(starts: number[], ends: number[], s: number, e: number) {
  let i = starts.length;
  while (i > 0 && starts[i - 1] > s) i--;
  starts.splice(i, 0, s);
  ends.splice(i, 0, e);
}

export class StickyLanes {
  private kept = new Map<Key, number>();
  private readonly starts: number[][] = [];
  private readonly ends: number[][] = [];

  /** The lanes kept, by item key (for tests and the viewer's debugging). */
  get size(): number {
    return this.kept.size;
  }

  laneOf(key: Key): number | undefined {
    return this.kept.get(key);
  }

  /**
   * Assigns a lane to each item (`order` lists the items by start), into `out`: −1 for one that does not fit in
   * `maxLanes`. With `settled`, sparse lanes are repacked. Returns the lanes used.
   */
  assign(
    keys: readonly Key[],
    starts: ArrayLike<number>,
    ends: ArrayLike<number>,
    order: ArrayLike<number>,
    maxLanes: number,
    settled: boolean,
    out: Int32Array,
  ): number {
    const sticky = this.pack(keys, starts, ends, order, maxLanes, out, true);
    if (settled && this.kept.size > 0) {
      const fresh = new Int32Array(out.length);
      const packed = this.pack(keys, starts, ends, order, maxLanes, fresh, false);
      if (packed.used < sticky.used && packed.overflow <= sticky.overflow) {
        out.set(fresh);
        return this.keep(keys, order, out, packed.used);
      }
    }
    return this.keep(keys, order, out, sticky.used);
  }

  private keep(keys: readonly Key[], order: ArrayLike<number>, lanes: Int32Array, used: number): number {
    const next = new Map<Key, number>();
    for (let k = 0; k < order.length; k++) {
      const i = order[k];
      if (lanes[i] >= 0) next.set(keys[i], lanes[i]);
    }
    this.kept = next;
    return used;
  }

  private pack(
    keys: readonly Key[],
    starts: ArrayLike<number>,
    ends: ArrayLike<number>,
    order: ArrayLike<number>,
    maxLanes: number,
    out: Int32Array,
    sticky: boolean,
  ): { used: number; overflow: number } {
    for (let l = 0; l < maxLanes; l++) {
      if (!this.starts[l]) {
        this.starts[l] = [];
        this.ends[l] = [];
      }
      this.starts[l].length = 0;
      this.ends[l].length = 0;
    }
    const extent = (i: number): [number, number] => [starts[i], Math.max(ends[i], starts[i] + 1)];
    let used = 0;
    out.fill(-1);
    if (sticky)
      for (let k = 0; k < order.length; k++) {
        const i = order[k];
        const l = this.kept.get(keys[i]);
        if (l === undefined || l >= maxLanes) continue;
        const [s, e] = extent(i);
        if (overlaps(this.starts[l], this.ends[l], s, e)) continue;
        insert(this.starts[l], this.ends[l], s, e);
        out[i] = l;
        if (l + 1 > used) used = l + 1;
      }
    let overflow = 0;
    for (let k = 0; k < order.length; k++) {
      const i = order[k];
      if (out[i] >= 0) continue;
      const [s, e] = extent(i);
      let lane = -1;
      for (let l = 0; l < maxLanes && lane < 0; l++) if (!overlaps(this.starts[l], this.ends[l], s, e)) lane = l;
      if (lane < 0) {
        overflow++;
        continue;
      }
      insert(this.starts[lane], this.ends[lane], s, e);
      out[i] = lane;
      if (lane + 1 > used) used = lane + 1;
    }
    return { used, overflow };
  }
}

/**
 * The hit index of marks in lanes: the marks of each lane sorted by their drawn left edge, with a running maximum of
 * right edges, so the marks under x are found in O(log n + what is there).
 */
export class LaneIndex {
  /** Mark indices by (lane, x0). */
  order = new Int32Array(0);
  /** Where each lane's marks start in `order`; `laneStart[lanes]` is the count. */
  laneStart = new Int32Array(1);
  runMax = new Float32Array(0);
  private x0: ArrayLike<number> = [];
  private x1: ArrayLike<number> = [];

  build(n: number, lanes: ArrayLike<number>, x0: ArrayLike<number>, x1: ArrayLike<number>, laneCount: number): void {
    this.x0 = x0;
    this.x1 = x1;
    const counts = new Int32Array(laneCount + 1);
    for (let i = 0; i < n; i++) if (lanes[i] >= 0) counts[lanes[i] + 1]++;
    for (let l = 0; l < laneCount; l++) counts[l + 1] += counts[l];
    this.laneStart = counts.slice();
    const fill = counts.slice();
    const order = new Int32Array(counts[laneCount]);
    for (let i = 0; i < n; i++) if (lanes[i] >= 0) order[fill[lanes[i]]++] = i;
    for (let l = 0; l < laneCount; l++) {
      const seg = order.subarray(counts[l], counts[l + 1]);
      seg.sort((a, b) => x0[a] - x0[b]);
    }
    const runMax = new Float32Array(order.length);
    for (let l = 0; l < laneCount; l++) {
      let m = Number.NEGATIVE_INFINITY;
      for (let k = counts[l]; k < counts[l + 1]; k++) {
        m = Math.max(m, x1[order[k]]);
        runMax[k] = m;
      }
    }
    this.order = order;
    this.runMax = runMax;
  }

  /** The mark under x in a lane (the last drawn there), or −1. */
  at(lane: number, x: number): number {
    if (lane < 0 || lane + 1 >= this.laneStart.length) return -1;
    const from = this.laneStart[lane];
    let lo = from;
    let hi = this.laneStart[lane + 1];
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.x0[this.order[mid]] <= x) lo = mid + 1;
      else hi = mid;
    }
    for (let k = lo - 1; k >= from && this.runMax[k] >= x; k--) {
      const i = this.order[k];
      if (this.x1[i] >= x) return i;
    }
    return -1;
  }
}
