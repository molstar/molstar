/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { TabixIndex } from '../tabix.js';

// A line's sequence and extent by a tabix index's rule (`_spec/bgzf.md` §3; htslib's `tbx_parse1`). Not published
// (`./common/internal/*`).

export const MaxPosition = 0x7fffffff;

/** An integer at `[s, e)` as `strtoll` reads it: optional sign, digits; NaN when none. */
export function parseInteger(bytes: Uint8Array, s: number, e: number): number {
  let i = s;
  let sign = 1;
  if (i < e && (bytes[i] === 43 || bytes[i] === 45)) {
    if (bytes[i] === 45) sign = -1;
    i++;
  }
  const d0 = i;
  let v = 0;
  while (i < e && bytes[i] >= 48 && bytes[i] <= 57 && v <= Number.MAX_SAFE_INTEGER / 10) v = v * 10 + bytes[i++] - 48;
  return i === d0 ? Number.NaN : sign * v;
}

/** The index of `pattern` in `bytes[s, e)`, or −1. */
function find(bytes: Uint8Array, s: number, e: number, pattern: readonly number[]): number {
  outer: for (let i = s; i + pattern.length <= e; i++) {
    for (let k = 0; k < pattern.length; k++) if (bytes[i + k] !== pattern[k]) continue outer;
    return i;
  }
  return -1;
}

const EndEq = [69, 78, 68, 61]; // END=
const SemiEndEq = [59, 69, 78, 68, 61]; // ;END=

export interface Extent {
  seqStart: number;
  seqEnd: number;
  begin: number;
  end: number;
}

/**
 * A line's sequence and extent by the index's rule (htslib's `tbx_parse1`): generic `begin` (−1 unless 0-based) to
 * `end`, or `begin + 1` without an end column; SAM to `POS − 1` plus the CIGAR's reference length; VCF to `INFO/END`
 * when INFO has `END=` past `POS − 1`, else `POS − 1 + len(REF)`. False when the line has no sequence or position.
 */
export function lineExtent(bytes: Uint8Array, s: number, e: number, index: TabixIndex, out: Extent): boolean {
  const { seq: cs, begin: cb, end: ce } = index.columns;
  const ucsc = index.zeroBased;
  let beg = -1;
  let end = -1;
  let seqStart = -1;
  let seqEnd = -1;
  for (let id = 1, b = s; b <= e; id++) {
    let i = bytes.indexOf(9, b);
    if (i < 0 || i > e) i = e;
    if (id === cs) {
      seqStart = b;
      seqEnd = i;
    } else if (id === cb) {
      const v = parseInteger(bytes, b, i);
      if (Number.isNaN(v)) return false;
      beg = v;
      if (cb <= ce) end = beg;
      if (!ucsc) beg -= 1;
      else if (cb <= ce) end += 1;
      if (beg < 0) beg = 0;
      if (end < 1) end = 1;
    } else if (index.format === 'generic') {
      if (id === ce) {
        const v = parseInteger(bytes, b, i);
        if (Number.isNaN(v)) return false;
        end = v;
      }
    } else if (index.format === 'sam') {
      if (id === 6 && beg >= 0) {
        let l = 0;
        for (let k = b; k < i; ) {
          let x = 0;
          while (k < i && bytes[k] >= 48 && bytes[k] <= 57) x = x * 10 + bytes[k++] - 48;
          const op = bytes[k++] | 0x20;
          if (op === 109 || op === 100 || op === 110 || op === 61 || op === 120) l += x; // M D N = X
        }
        end = beg + (l === 0 ? 1 : l);
      }
    } else if (id === 4) {
      if (b < i && beg >= 0) end = beg + (i - b);
    } else if (id === 8 && beg >= 0) {
      let at = -1;
      if (find(bytes, b, i, EndEq) === b) at = b + 4;
      else {
        const semi = find(bytes, b, i, SemiEndEq);
        if (semi >= 0) at = semi + 5;
      }
      if (at >= 0 && at < i && bytes[at] !== 46) {
        const v = parseInteger(bytes, at, i);
        if (v > beg) end = v;
      }
    }
    if (i >= e) break;
    b = i + 1;
  }
  if (seqStart < 0 || beg < 0) return false;
  // without an end column (`end` 0) a generic line covers its position
  if (end < 0 || (index.format !== 'vcf' && index.format !== 'sam' && ce === 0)) end = beg + 1;
  out.seqStart = seqStart;
  out.seqEnd = seqEnd;
  out.begin = beg;
  out.end = end;
  return true;
}

export function sameBytes(bytes: Uint8Array, s: number, e: number, name: Uint8Array) {
  if (e - s !== name.length) return false;
  for (let k = 0; k < name.length; k++) if (bytes[s + k] !== name[k]) return false;
  return true;
}
