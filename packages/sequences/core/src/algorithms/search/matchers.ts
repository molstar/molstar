/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The sequence matchers of `@molstar/seq-core/algorithms/search` (seqspec:libraries#search): the contract every matcher shares,
// the exact matcher (Horspool), the class matcher (shift-and) and the motif compiler. Hits go to a callback; nothing is
// allocated per hit, and no matcher knows an alphabet: a code is a byte, and what it stands for is the caller's.

/** Called once per hit: `[start, stop)` in the searched array's indices (`start > stop` when it crosses the end of a cyclic search), `score` the edit count (0 but for `approximateMatcher`). */
export type HitVisitor = (start: number, stop: number, score: number) => void;

export interface Matcher {
  /**
   * The longest hit the matcher can report: the pattern's length, the pattern's plus `maxEdits` for the approximate
   * matcher, `maxLength` for the regex one. A cyclic search reads this many codes past `to`, so a hit may cross the end.
   */
  readonly length: number;
  /**
   * Calls `visit` for every hit starting in `[from, to)` of `seq`, in order of their ends, and returns how many. A
   * hit lies inside `[from, to)` (it ends by `to`); with `cyclic`, `seq` is a circle and a hit may run on past `to`,
   * across the end into the codes from 0 (`start > stop`, as the format writes a wrapping range,
   * seqspec:format#positions; a hit as long as the circle has `start == stop`, unless it starts at 0: `[0,
   * seq.length)`). An empty pattern has no hits. Hits may overlap. Not reentrant: `visit` must not call `find` of the
   * same matcher. Throws a `RangeError` for a window outside `[0, seq.length]`.
   */
  find(seq: Uint8Array, from: number, to: number, cyclic: boolean, visit: HitVisitor): number;
}

/** Throws a `RangeError` unless `0 <= from <= to <= length`, all integers. */
export function checkWindow(length: number, from: number, to: number): void {
  if (!Number.isInteger(from) || !Number.isInteger(to) || from < 0 || to > length || from > to) {
    throw new RangeError(`search: the window [${from}, ${to}) must be integers with 0 <= from <= to <= ${length}`);
  }
}

/**
 * One past the last code a search of `[from, to)` reads. Not cyclic: `to`. Cyclic: the hits' starts end at `to`, so the
 * last hit of `matcherLength` codes ends `matcherLength - 1` later, but a hit is at most one lap long. Positions at or
 * beyond `seq.length` stand for `position - seq.length`.
 */
export function readEnd(matcherLength: number, seqLength: number, to: number, cyclic: boolean): number {
  return cyclic ? Math.min(to + matcherLength - 1, to + seqLength) : to;
}

/** The codes of a string's UTF-16 units, for strings of Latin-1 letters. Throws a `RangeError` for a unit above 255. */
export function toCodes(text: string): Uint8Array {
  const out = new Uint8Array(text.length);
  for (let i = 0; i < text.length; i++) {
    const c = text.charCodeAt(i);
    if (c > 255)
      throw new RangeError(`search: "${text[i]}" (U+${c.toString(16).toUpperCase()}) is not a code from 0 to 255`);
    out[i] = c;
  }
  return out;
}

/** Horspool's search for one exact pattern. The pattern is copied. */
export function exactMatcher(pattern: Uint8Array): Matcher {
  const pat = pattern.slice();
  const m = pat.length;
  const shift = new Int32Array(256).fill(m);
  for (let i = 0; i < m - 1; i++) shift[pat[i]] = m - 1 - i;
  return {
    length: m,
    find(seq, from, to, cyclic, visit) {
      checkWindow(seq.length, from, to);
      const n = seq.length;
      // A pattern longer than the circle is not a hit on it.
      if (m === 0 || from >= to || (cyclic && m > n)) return 0;
      const limit = readEnd(m, n, to, cyclic);
      const last = pat[m - 1];
      let count = 0;
      // The hits inside the array: Horspool.
      const inside = Math.min(to, Math.min(limit, n) - m + 1);
      for (let s = from; s < inside; ) {
        const c = seq[s + m - 1];
        if (c === last) {
          let k = m - 2;
          while (k >= 0 && seq[s + k] === pat[k]) k--;
          if (k < 0) {
            visit(s, s + m, 0);
            count++;
          }
        }
        s += shift[c];
      }
      // The hits across the end of a circle, checked at each start from the end.
      if (cyclic) {
        for (let s = Math.max(from, n - m + 1); s < to && s + m <= limit; s++) {
          let k = 0;
          while (k < m && seq[s + k >= n ? s + k - n : s + k] === pat[k]) k++;
          if (k === m) {
            visit(s, s + m - n, 0);
            count++;
          }
        }
      }
      return count;
    },
  };
}

/** The most positions of a `classMatcher` pattern (two 32-bit words). */
export const MAX_CLASS_PATTERN = 64;

/**
 * Shift-and for a pattern of sets: position i matches any code of `pattern[i]`. Patterns of up to 64 positions
 * (`MAX_CLASS_PATTERN`): longer ones throw a `RangeError`; use `approximateMatcher` or `regexMatcher` for them.
 */
export function classMatcher(pattern: readonly Uint8Array[]): Matcher {
  const m = pattern.length;
  if (m > MAX_CLASS_PATTERN) {
    throw new RangeError(
      `classMatcher: ${m} positions; the limit is ${MAX_CLASS_PATTERN} (use regexMatcher for longer patterns)`,
    );
  }
  // maskLo/maskHi[c]: bit i set when position i (the low word: 0 to 31; the high word: 32 to 63) accepts code c.
  const maskLo = new Int32Array(256);
  const maskHi = new Int32Array(256);
  for (let i = 0; i < m; i++) {
    for (const c of pattern[i]) {
      if (i < 32) maskLo[c] |= 1 << i;
      else maskHi[c] |= 1 << (i - 32);
    }
  }
  const acceptLo = m >= 1 && m <= 32 ? 1 << (m - 1) : 0;
  const acceptHi = m > 32 ? 1 << (m - 33) : 0;
  return {
    length: m,
    find(seq, from, to, cyclic, visit) {
      checkWindow(seq.length, from, to);
      const n = seq.length;
      if (m === 0 || from >= to || (cyclic && m > n)) return 0;
      const limit = readEnd(m, n, to, cyclic);
      let lo = 0;
      let hi = 0;
      let count = 0;
      for (let u = from; u < limit; u++) {
        const c = seq[u >= n ? u - n : u];
        hi = ((hi << 1) | (lo >>> 31)) & maskHi[c];
        lo = ((lo << 1) | 1) & maskLo[c];
        if (((lo & acceptLo) | (hi & acceptHi)) !== 0) {
          visit(u + 1 - m, u + 1 > n ? u + 1 - n : u + 1, 0);
          count++;
        }
      }
      return count;
    },
  };
}

/**
 * A pattern of code sets from letters: each letter of `motif` stands for the letters `classes` gives it (`R: "AG"`), or
 * for itself when `classes` has no entry. Letters are Latin-1 (codes 0 to 255), case-sensitive, and the sets are the
 * caller's data: IUPAC's nucleotide codes, a grouping of amino acids. Throws a `RangeError` for a letter above 255.
 */
export function compileMotif(motif: string, classes: Readonly<Record<string, string>>): readonly Uint8Array[] {
  const out: Uint8Array[] = [];
  for (const letter of motif) {
    const set = Object.hasOwn(classes, letter) ? classes[letter] : letter;
    out.push(Uint8Array.from(new Set(toCodes(set))).sort());
  }
  return out;
}
