/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The text index of `@molstar/seq-core/algorithms/search` (seqspec:libraries#search): sorted normalized keys with their ids, for
// the names, labels and fields a provider streams.

export interface TextIndex {
  /** Adds `texts[start + i]` for the id `ids[start + i]`, for `i < count`. An id may have several texts; a text added again for its id is held once, so repeating as data streams in is safe. */
  add(ids: Int32Array, texts: readonly string[], start: number, count: number): void;
  /** Removes every text of the first `count` ids of `ids`. */
  remove(ids: Int32Array, count: number): void;
  /** Visits the ids of the texts equal to `q` after normalization, in key order (an id once per distinct matching text); returns the count. */
  exact(q: string, visit: (id: number) => void): number;
  /** The same for the texts that start with `q` (every text for an empty `q`): two binary searches. */
  prefix(q: string, visit: (id: number) => void): number;
  /** The same for the normalized texts `re` matches anywhere: a scan. The expression is written for normalized text (lower case by default). */
  matches(re: RegExp, visit: (id: number) => void): number;
  /** The number of distinct (id, text) pairs held. */
  readonly size: number;
}

export interface TextIndexOptions {
  /** What keys and queries pass through; default NFKC and lower case. */
  readonly normalize?: (s: string) => string;
}

const byKey = (a: string, b: string): number => (a < b ? -1 : a > b ? 1 : 0);
// NFKC and `toLowerCase`: not full case folding ("STRASSE" does not find "Straße"); a host that needs it passes `normalize`.
const defaultNormalize = (s: string): string => s.normalize("NFKC").toLowerCase();

class SortedTextIndex implements TextIndex {
  private readonly normalize: (s: string) => string;
  /** The sorted keys with their ids (ties by id), and the keys added since, merged in when a query needs them. */
  private keys: string[] = [];
  private ids = new Int32Array(0);
  private pendingKeys: string[] = [];
  private pendingIds: number[] = [];

  constructor(normalize: (s: string) => string) {
    this.normalize = normalize;
  }

  get size(): number {
    this.merge();
    return this.keys.length;
  }

  add(ids: Int32Array, texts: readonly string[], start: number, count: number): void {
    for (let i = start; i < start + count; i++) {
      this.pendingKeys.push(this.normalize(texts[i]));
      this.pendingIds.push(ids[i]);
    }
  }

  remove(ids: Int32Array, count: number): void {
    if (count === 0 || this.size === 0) return;
    const gone = new Set<number>();
    for (let i = 0; i < count; i++) gone.add(ids[i]);
    this.merge();
    const keys: string[] = [];
    const kept: number[] = [];
    for (let i = 0; i < this.keys.length; i++) {
      if (gone.has(this.ids[i])) continue;
      keys.push(this.keys[i]);
      kept.push(this.ids[i]);
    }
    this.keys = keys;
    this.ids = Int32Array.from(kept);
  }

  exact(q: string, visit: (id: number) => void): number {
    this.merge();
    const key = this.normalize(q);
    let count = 0;
    for (let i = this.lowerBound(key); i < this.keys.length && this.keys[i] === key; i++, count++) visit(this.ids[i]);
    return count;
  }

  prefix(q: string, visit: (id: number) => void): number {
    this.merge();
    const key = this.normalize(q);
    // The keys from the first one not below `key` that start with it come first and in a row.
    const first = this.lowerBound(key);
    let lo = first;
    let hi = this.keys.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (this.keys[mid].startsWith(key)) lo = mid + 1;
      else hi = mid;
    }
    for (let i = first; i < lo; i++) visit(this.ids[i]);
    return lo - first;
  }

  matches(re: RegExp, visit: (id: number) => void): number {
    this.merge();
    const stateful = re.global || re.sticky;
    let count = 0;
    for (let i = 0; i < this.keys.length; i++) {
      if (stateful) re.lastIndex = 0;
      if (!re.test(this.keys[i])) continue;
      visit(this.ids[i]);
      count++;
    }
    return count;
  }

  /** The index of the first key not below `key`. */
  private lowerBound(key: string): number {
    let lo = 0;
    let hi = this.keys.length;
    while (lo < hi) {
      const mid = (lo + hi) >>> 1;
      if (this.keys[mid] < key) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  /** Sorts the keys added since the last query and merges them into the sorted ones. */
  private merge(): void {
    const p = this.pendingKeys.length;
    if (p === 0) return;
    const order = Array.from({ length: p }, (_, i) => i).sort(
      (a, b) => byKey(this.pendingKeys[a], this.pendingKeys[b]) || this.pendingIds[a] - this.pendingIds[b],
    );
    const total = this.keys.length + p;
    // Pre-sized with a write cursor, truncated at the end (LA-08: pushing a tuple per element cost 4× the time).
    const keys: string[] = new Array(total);
    const ids = new Int32Array(total);
    let w = 0;
    let a = 0;
    let b = 0;
    for (let o = 0; o < total; o++) {
      const takeOld =
        b >= p ||
        (a < this.keys.length &&
          (byKey(this.keys[a], this.pendingKeys[order[b]]) || this.ids[a] - this.pendingIds[order[b]]) <= 0);
      let key: string;
      let id: number;
      if (takeOld) {
        key = this.keys[a];
        id = this.ids[a++];
      } else {
        key = this.pendingKeys[order[b]];
        id = this.pendingIds[order[b++]];
      }
      // GM-08: the same text of the same id, added again (a region fetched again), is held once.
      if (w > 0 && keys[w - 1] === key && ids[w - 1] === id) continue;
      keys[w] = key;
      ids[w++] = id;
    }
    keys.length = w;
    this.keys = keys;
    this.ids = w === total ? ids : ids.slice(0, w);
    this.pendingKeys = [];
    this.pendingIds = [];
  }
}

/** A text index; `normalize` defaults to NFKC and lower case. */
export function textIndex(opts?: TextIndexOptions): TextIndex {
  return new SortedTextIndex(opts?.normalize ?? defaultNormalize);
}
