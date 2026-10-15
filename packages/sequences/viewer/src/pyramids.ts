/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The index for scale (seqspec:sources#levels-and-the-zoom-gate): segment tables and summary pyramids, so a frame
// reduces O(device columns) bins and visits only the items it draws, at any zoom. Segments: start, stop sorted by start
// with a prefix maximum of the stops, a window being two predecessor searches (the legacy Seq* viewer's segment table,
// after Mol*'s `findPredecessorIndexInRange`, written fresh). Pyramids: per block, power-of-two bins over model indices
// from `PYRAMID_BASE` up; below it a draw reads raw positions (at most 2^BASE per device column). Built once (prepare,
// or a draw's first use per theme for classes and counts, which are colors), never per frame; queries allocate nothing.

/** The last index i in [lo, hi) with a[i] <= x; lo - 1: none. a is non-decreasing there. */
export function predecessor(a: Int32Array, lo: number, hi: number, x: number): number {
  let l = lo;
  let h = hi;
  while (l < h) {
    const m = (l + h) >> 1;
    if (a[m] <= x) l = m + 1;
    else h = m;
  }
  return l - 1;
}

/**
 * Per segment [lo, hi) of items sorted by start: `reach[k]`, the largest stop of items lo..k (a boundary one position
 * wide), so items reaching `from` start at the first k whose reach passes it, overlapping or not.
 */
export function prefixReach(start: Int32Array, stop: Int32Array, lo: number, hi: number, out: Int32Array) {
  let r = -2147483648;
  for (let k = lo; k < hi; k++) {
    r = Math.max(r, stop[k], start[k] + 1);
    out[k] = r;
  }
}
/**
 * The items of segment [lo, hi) that may overlap [from, to) into `out` ([out[0], out[1])): those from the first whose
 * reach passes `from` to the first starting at `to` or later. Overlapping items in it may still end before `from`.
 */
export function segWindow(
  start: Int32Array,
  reach: Int32Array,
  lo: number,
  hi: number,
  from: number,
  to: number,
  out: Int32Array,
) {
  out[0] = predecessor(reach, lo, hi, from) + 1;
  out[1] = Math.max(out[0], predecessor(start, lo, hi, to - 1) + 1);
}

/** The finest level kept: bins of 2^BASE positions. */
export const PYRAMID_BASE = 2;
/**
 * The level (log2 of its bin) for `per` positions per device column: the largest bin not wider, of `levels` from
 * `base`; -1 below the base.
 */
export function levelFor(per: number, levels: number, base = PYRAMID_BASE): number {
  if (!(per >= 2 ** base) || levels === 0) return -1;
  return Math.min(Math.floor(Math.log2(per)), base + levels - 1);
}

/**
 * Per-position classes (1–255, 0: none) by level: per bin the plurality class (0 only where nothing else) and its share
 * of the bin × 255. The base is exact; above it a bin takes the larger of its halves' classes (their counts summed
 * when they agree), as a tournament: a class holding more than half a bin is found exactly, otherwise the class is the
 * halves' plurality; the share is a lower bound (exact where both halves agree all the way down). Ties go by a hash of the bin, not to the lower
 * class: a tournament amplifies any fixed preference (a nucleotide track read as all A).
 */
export interface ClassPyramid {
  readonly cls: readonly Uint8Array[];
  readonly share: readonly Uint8Array[];
}
const tally = new Uint32Array(256);
/** A bin's tie-break: 0 or 1, by a multiplicative hash of its index and level. */
const coin = (j: number, k: number) => (Math.imul(j + 1, 0x9e3779b1) ^ Math.imul(k + 1, 0x85ebca6b)) >>> 31;
export function classPyramid(classes: Uint8Array): ClassPyramid {
  return classPyramidFrom(classes.length, (s, out) => {
    out.set(classes.subarray(s, s + out.length));
    return true;
  }) as ClassPyramid;
}
/** The base level is filled 64 KiB at a time (never an array a block long beside the pyramid). */
const CHUNK = 1 << 16;
const chunk = new Uint8Array(CHUNK);
/**
 * `classPyramid` over n positions whose classes `fill(start, out)` writes, `out.length` at a time from `start`; a fill
 * that returns false stops the build (null): more classes than a byte holds.
 */
export function classPyramidFrom(n: number, fill: (start: number, out: Uint8Array) => boolean): ClassPyramid | null {
  const B = 1 << PYRAMID_BASE;
  const cls: Uint8Array[] = [];
  const share: Uint8Array[] = [];
  let bins = Math.ceil(n / B);
  let c = new Uint8Array(bins);
  let cnt = new Uint32Array(bins);
  for (let s = 0; s < n; s += CHUNK) {
    const part = chunk.subarray(0, Math.min(CHUNK, n - s));
    if (!fill(s, part)) return null;
    for (let j = s / B; j * B < s + part.length; j++) {
      const e = Math.min(part.length, (j + 1) * B - s);
      let best = 0;
      for (let i = j * B - s; i < e; i++) {
        const x = part[i];
        if (x === 0) continue;
        tally[x]++;
        if (best === 0 || tally[x] > tally[best] || (tally[x] === tally[best] && coin(j, i + s) === 1)) best = x;
      }
      c[j] = best;
      cnt[j] = best > 0 ? tally[best] : 0;
      for (let i = j * B - s; i < e; i++) tally[part[i]] = 0;
    }
  }
  for (let size = B, k = PYRAMID_BASE; ; size *= 2, k++) {
    cls.push(c);
    const sh = new Uint8Array(bins);
    for (let j = 0; j < bins; j++) sh[j] = Math.round((cnt[j] / size) * 255);
    share.push(sh);
    if (bins <= 1) break;
    const up = Math.ceil(bins / 2);
    const c2 = new Uint8Array(up);
    const n2 = new Uint32Array(up);
    for (let j = 0; j < up; j++) {
      const a = 2 * j;
      const b = a + 1 < bins ? a + 1 : a;
      const pick =
        b !== a && (cnt[b] > cnt[a] || (cnt[b] === cnt[a] && c[b] !== 0 && (c[a] === 0 || coin(j, k) === 1)));
      c2[j] = pick ? c[b] : c[a];
      n2[j] = b !== a && c[a] === c[b] ? cnt[a] + cnt[b] : pick ? cnt[b] : cnt[a];
    }
    c = c2;
    cnt = n2;
    bins = up;
  }
  return { cls, share };
}

/**
 * Values by level: per bin the minimum, maximum, sum and count of the values that are not NaN (none: NaN, 0). Doubles,
 * as the values (seqspec:format#values): a bin of one value holds that value exactly. Without `minmax` (a `mean`
 * theme's reduction) only the sums and counts: `min` and `max` are empty.
 */
export interface ValuePyramid {
  readonly min: readonly Float64Array[];
  readonly max: readonly Float64Array[];
  readonly sum: readonly Float64Array[];
  readonly count: readonly Uint32Array[];
}
export function valuePyramid(values: Float64Array, minmax = true): ValuePyramid {
  const n = values.length;
  const B = 1 << PYRAMID_BASE;
  const out = {
    min: [] as Float64Array[],
    max: [] as Float64Array[],
    sum: [] as Float64Array[],
    count: [] as Uint32Array[],
  };
  let bins = Math.ceil(n / B);
  let lo = new Float64Array(minmax ? bins : 0).fill(Number.NaN);
  let hi = new Float64Array(minmax ? bins : 0).fill(Number.NaN);
  let sum = new Float64Array(bins);
  let cnt = new Uint32Array(bins);
  for (let i = 0; i < n; i++) {
    const x = values[i];
    if (Number.isNaN(x)) continue;
    const j = i >> PYRAMID_BASE;
    if (minmax && !(x >= lo[j])) lo[j] = x;
    if (minmax && !(x <= hi[j])) hi[j] = x;
    sum[j] += x;
    cnt[j]++;
  }
  for (;;) {
    if (minmax) {
      out.min.push(lo);
      out.max.push(hi);
    }
    out.sum.push(sum);
    out.count.push(cnt);
    if (bins <= 1) return out;
    const up = Math.ceil(bins / 2);
    const m = minmax ? up : 0;
    const [l2, h2, s2, c2] = [new Float64Array(m), new Float64Array(m), new Float64Array(up), new Uint32Array(up)];
    for (let j = 0; j < up; j++) {
      const a = 2 * j;
      const b = a + 1 < bins ? a + 1 : a;
      if (minmax) {
        l2[j] = cnt[a] === 0 ? lo[b] : cnt[b] === 0 ? lo[a] : Math.min(lo[a], lo[b]);
        h2[j] = cnt[a] === 0 ? hi[b] : cnt[b] === 0 ? hi[a] : Math.max(hi[a], hi[b]);
      }
      s2[j] = sum[a] + (b !== a ? sum[b] : 0);
      c2[j] = cnt[a] + (b !== a ? cnt[b] : 0);
    }
    [lo, hi, sum, cnt, bins] = [l2, h2, s2, c2, up];
  }
}

/**
 * Items by level, from `base`: per bin and slot (an item's color class, < `slots`) the items overlapping the bin, or
 * with `point` each item once, in the bin of its center (points and variants: a codon-long point counts once),
 * bin-major (`counts[level][bin * slots + slot]`). Positions are `start - offset` (a block's domain start). Each level
 * from a difference array: O(items + bins). Float32: drawing counts, never data values (exact below 2^24 items a bin).
 */
export interface CountPyramid {
  readonly base: number;
  readonly slots: number;
  readonly counts: readonly Float32Array[];
}
export interface CountItems {
  readonly start: Int32Array;
  readonly stop: Int32Array;
  readonly slot: Uint8Array;
  readonly from: number;
  readonly to: number;
  readonly offset?: number;
  readonly point?: boolean;
}
export function countPyramid(length: number, base: number, slots: number, items: CountItems): CountPyramid {
  const counts: Float32Array[] = [];
  const { start, stop, slot, from, to } = items;
  const o = items.offset ?? 0;
  for (let k = base; ; k++) {
    const bins = Math.ceil(length / 2 ** k);
    const d = new Float32Array((bins + 1) * slots);
    for (let i = from; i < to; i++) {
      const e = Math.max(stop[i], start[i] + 1) - 1 - o;
      // Unsigned: a sum past 2^31 (a block of over 1.07e9 positions) stays positive (VL-18).
      const c = (start[i] - o + e) >>> 1;
      const j0 = Math.max(0, Math.floor((items.point ? c : start[i] - o) / 2 ** k));
      const j1 = Math.min(bins - 1, Math.floor((items.point ? c : e) / 2 ** k));
      if (j1 < j0) continue;
      d[j0 * slots + slot[i]]++;
      d[(j1 + 1) * slots + slot[i]]--;
    }
    for (let j = 1; j <= bins; j++) for (let s = 0; s < slots; s++) d[j * slots + s] += d[(j - 1) * slots + s];
    counts.push(d.subarray(0, bins * slots));
    if (bins <= 1) return { base, slots, counts };
  }
}

// A draw's reads: a device column takes the bins of a level whose start lies in it, [j0, j1), so a copied strip and a
// full redraw count the same bins. Allocation-free: results into `out`, never a double returned.

const acc = new Float64Array(256);
const seen = new Uint8Array(256);
const mark = new Uint8Array(256);
/**
 * The plurality class of bins [j0, j1) at `level` by their summed shares (ties by a hash of the column, as the
 * tournament's), 0: none; `out[0]` its summed share over the bins' positions, 0..255.
 */
export function binClass(p: ClassPyramid, level: number, j0: number, j1: number, out: Float64Array): number {
  const L = level - PYRAMID_BASE;
  const cls = p.cls[L];
  const sh = p.share[L];
  const hi = Math.min(j1, cls.length);
  let n = 0;
  let best = 0;
  for (let j = Math.max(0, j0); j < hi; j++) {
    const c = cls[j];
    if (c === 0) continue;
    if (mark[c] === 0) {
      mark[c] = 1;
      seen[n++] = c;
    }
    acc[c] += sh[j];
    if (best === 0 || acc[c] > acc[best] || (acc[c] === acc[best] && c !== best && coin(j0, c) === 1)) best = c;
  }
  out[0] = best > 0 && hi > j0 ? acc[best] / (hi - Math.max(0, j0)) : 0;
  for (let k = 0; k < n; k++) acc[seen[k]] = mark[seen[k]] = 0;
  return best;
}
/** Bins [j0, j1) at `level`: `out` the minimum, maximum (NaN without them), sum and count of their values. */
export function binValues(p: ValuePyramid, level: number, j0: number, j1: number, out: Float64Array) {
  const L = level - PYRAMID_BASE;
  const sum = p.sum[L];
  const cnt = p.count[L];
  const mm = p.min.length > 0;
  out[0] = out[1] = Number.NaN;
  out[2] = out[3] = 0;
  for (let j = Math.max(0, j0), hi = Math.min(j1, sum.length); j < hi; j++) {
    if (cnt[j] === 0) continue;
    if (mm && !(p.min[L][j] >= out[0])) out[0] = p.min[L][j];
    if (mm && !(p.max[L][j] <= out[1])) out[1] = p.max[L][j];
    out[2] += sum[j];
    out[3] += cnt[j];
  }
}
/** Bins [j0, j1) at `level`: per slot into `out` the sum of their counts, or the maximum (`max`: coverage). */
export function binCounts(p: CountPyramid, level: number, j0: number, j1: number, max: boolean, out: Float64Array) {
  const c = p.counts[level - p.base];
  const S = p.slots;
  out.fill(0, 0, S);
  for (let j = Math.max(0, j0), hi = Math.min(j1, c.length / S); j < hi; j++)
    for (let s = 0; s < S; s++) {
      const x = c[j * S + s];
      out[s] = max ? (x > out[s] ? x : out[s]) : out[s] + x;
    }
}
