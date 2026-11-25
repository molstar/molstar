/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The binning scheme of TBI and CSI (SAM spec §5.3, CSIv1 spec; `_spec/bgzf.md`): bins of a region, a bin's level
// and parent. Not published (`./common/internal/*`); tests reach it here. Arithmetic is by division, not shifts: a CSI
// may address past 2^31 (`minShift + 3 × depth` up to 52).

/** The first bin of level `level` (0 the whole range): (8^level − 1) / 7. */
export function binFirst(level: number): number {
  return (8 ** level - 1) / 7;
}

/** The bins of a scheme with `depth` levels below the root: bins `0 … binFirst(depth + 1) − 1`. */
export function binCount(depth: number): number {
  return binFirst(depth + 1);
}

/** The pseudo-bin of an index's metadata (TBI: 37,450; its bins end at 37,448): `binCount + 1`, as htslib's. */
export function metaBin(depth: number): number {
  return binCount(depth) + 1;
}

/** The parent of bin `bin` (> 0). */
export function binParent(bin: number): number {
  return Math.floor((bin - 1) / 8);
}

/** The smallest bin holding `[beg, end)` (0-based, half-open; `end > beg`): htslib's `hts_reg2bin`. */
export function reg2bin(beg: number, end: number, minShift: number, depth: number): number {
  end -= 1;
  let s = minShift;
  let t = binFirst(depth);
  for (let l = depth; l > 0; l--) {
    const size = 2 ** s;
    if (Math.floor(beg / size) === Math.floor(end / size)) return t + Math.floor(beg / size);
    s += 3;
    t -= 8 ** (l - 1);
  }
  return 0;
}

/** Every bin that may hold a record meeting `[beg, end)` (0-based, half-open), by level: htslib's `hts_reg2bins`. */
export function reg2bins(beg: number, end: number, minShift: number, depth: number): number[] {
  const bins: number[] = [];
  if (beg < 0) beg = 0;
  let s = minShift + 3 * depth;
  const max = 2 ** s;
  if (end > max) end = max;
  if (beg >= end) return bins;
  end -= 1;
  for (let l = 0, t = 0; l <= depth; l++) {
    const size = 2 ** s;
    const b = t + Math.floor(beg / size);
    const e = t + Math.floor(end / size);
    for (let i = b; i <= e; i++) bins.push(i);
    s -= 3;
    t += 8 ** l;
  }
  return bins;
}
