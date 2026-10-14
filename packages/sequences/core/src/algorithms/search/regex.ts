/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The regular expression matcher of `@molstar/seq-core/algorithms/search` (seqspec:libraries#search): a JavaScript RegExp over a
// Latin-1 string view of the codes.

import { checkWindow, type HitVisitor, type Matcher, readEnd } from "./matchers.js";

/** The most code units `String.fromCharCode` is spread over at once. */
const CHUNK = 8192;

/** The codes at the unrolled positions `[from, limit)` of `seq` (positions from `seq.length` on wrap) as one string. */
function latin1(seq: Uint8Array, from: number, limit: number): string {
  const n = seq.length;
  let out = "";
  for (let u = from; u < limit; ) {
    const stop = Math.min(limit, u < n ? n : Number.POSITIVE_INFINITY, u + CHUNK);
    const a = u >= n ? u - n : u;
    out += String.fromCharCode(...seq.subarray(a, a + (stop - u)));
    u = stop;
  }
  return out;
}

/**
 * Pending hits: a binary min-heap by (stop, start) of `size` over two arrays that only grow, a find's own (a visitor
 * may search again).
 */
interface Pending {
  readonly start: number[];
  readonly stop: number[];
  size: number;
}
function push(h: Pending, start: number, stop: number): void {
  const S = h.start;
  const T = h.stop;
  let k = h.size++;
  if (k === S.length) {
    S.push(start);
    T.push(stop);
  }
  while (k > 0) {
    const up = (k - 1) >> 1;
    if (T[up] < stop || (T[up] === stop && S[up] < start)) break;
    S[k] = S[up];
    T[k] = T[up];
    k = up;
  }
  S[k] = start;
  T[k] = stop;
}
/** Removes the least: the last moved to the root and sifted down. */
function pop(h: Pending): void {
  const S = h.start;
  const T = h.stop;
  const last = --h.size;
  const start = S[last];
  const stop = T[last];
  let k = 0;
  for (;;) {
    let m = 2 * k + 1;
    if (m >= last) break;
    if (m + 1 < last && (T[m + 1] < T[m] || (T[m + 1] === T[m] && S[m + 1] < S[m]))) m++;
    if (stop < T[m] || (stop === T[m] && start < S[m])) break;
    S[k] = S[m];
    T[k] = T[m];
    k = m;
  }
  S[k] = start;
  T[k] = stop;
}

/**
 * A regular expression over the codes read as Latin-1 letters. Every start is tried, so hits overlap (`AA` in `AAA` is
 * two) and the leftmost, JavaScript-preferred match from each start is the hit; empty matches are not hits, and nor
 * is a match longer than a cyclic search's circle. The
 * expression sees only the window (and, cyclic, `maxLength - 1` codes beyond it): `maxLength` is the longest match to find
 * across the end of a cyclic search, 64 by default, and what `length` reports. The one matcher that allocates per hit:
 * the engine's match array, and the overlapping hits held to be visited in order of their ends, as every matcher's
 * are (the engine finds them by start; CB-04). Builds a string of the window for each search.
 */
export function regexMatcher(re: RegExp, maxLength = 64): Matcher {
  if (!Number.isInteger(maxLength) || maxLength < 1)
    throw new RangeError(`regexMatcher: maxLength must be an integer of at least 1; got ${maxLength}`);
  const rx = new RegExp(re.source, re.flags.includes("g") ? re.flags : `${re.flags}g`);
  return {
    length: maxLength,
    find(seq: Uint8Array, from: number, to: number, cyclic: boolean, visit: HitVisitor): number {
      checkWindow(seq.length, from, to);
      if (from >= to) return 0;
      const n = seq.length;
      const text = latin1(seq, from, readEnd(maxLength, n, to, cyclic));
      // The hits pending by end (then start: one hit a start), a binary min-heap: the engine yields starts in
      // increasing order and a hit ends past its start, so before a match at `i` every pending hit ending by `i + 1`
      // is due (LA-09: memory bounded by the overlap, where sorting all of them held every hit).
      const h: Pending = { start: [], stop: [], size: 0 };
      let count = 0;
      rx.lastIndex = 0;
      for (let match = rx.exec(text); ; match = rx.exec(text)) {
        const i = match === null ? -1 : match.index;
        const done = match === null || i >= to - from;
        // Due: every pending hit ending by `from + i + 1` (all of them at the end).
        const limit = done ? Number.POSITIVE_INFINITY : from + i + 1;
        while (h.size > 0 && h.stop[0] <= limit) {
          const start = h.start[0];
          const stop = h.stop[0];
          pop(h);
          visit(start, stop > n ? stop - n : stop, 0);
          count++;
        }
        if (done) break;
        rx.lastIndex = i + 1;
        const length = (match as RegExpExecArray)[0].length;
        // Empty matches are not hits, and nor is a match longer than the circle it wraps.
        if (length === 0 || length > n) continue;
        push(h, from + i, from + i + length);
      }
      return count;
    },
  };
}
