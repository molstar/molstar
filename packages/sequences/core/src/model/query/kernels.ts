/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The selection builder and the run kernels (seqspec:seqql#the-selection, seqspec:seqql#performance). A builder grows
// its arrays by doubling and is reused within a run of a query. Runs within a group are in canonical order (system:
// kind, id, row; then start), disjoint, merged where they touch. The kernels work on runs, never per position, and
// allocate nothing.

import { CoordinateSystemKind as K } from "../coordinates.ts";
import type { SeqSet } from "../set.ts";
import type { SeqSelection } from "./selection.ts";

/** What a builder reports to: the context's limits (runs, groups) and its clock, read every 1,024 runs. */
export interface BuildGuard {
  readonly maxRuns: number;
  readonly maxGroups: number;
  over(what: "runs" | "groups", count: number): never;
  tick(): void;
}

type Runs = Pick<SeqSelection, "kind" | "id" | "row" | "start" | "stop">;
const grow = <T extends Int32Array | Uint8Array | Float64Array>(a: T, n: number): T => {
  const b = new (a.constructor as new (n: number) => T)(Math.max(n, 2 * a.length));
  b.set(a);
  return b;
};
/** Run j's system against (kind, id, row): negative, zero or positive. */
const cmpSys = (s: Runs, j: number, kind: number, id: number, row: number) =>
  s.kind[j] - kind || s.id[j] - id || s.row[j] - row;

export class SelectionBuilder {
  runCount = 0;
  groupCount = 0;
  kind = new Uint8Array(64);
  id = new Int32Array(64);
  row = new Int32Array(64);
  start = new Int32Array(64);
  stop = new Int32Array(64);
  groupStart = new Int32Array(17);
  groupAnnotation = new Int32Array(16);
  groupItem = new Int32Array(16);
  groupScore = new Float64Array(16);
  guard: BuildGuard | null = null;
  private open = -1;
  private sorted = true;

  reset(guard: BuildGuard | null): this {
    this.runCount = this.groupCount = 0;
    this.open = -1;
    this.guard = guard;
    return this;
  }

  begin(annotation = -1, item = -1, score = Number.NaN): void {
    const g = this.groupCount;
    // Grown by the per-group arrays, so `groupStart` (one longer) never outgrows them (QR2a's oracle: tags were lost).
    if (g >= this.groupAnnotation.length) {
      this.groupStart = grow(this.groupStart, 2 * g + 2);
      this.groupAnnotation = grow(this.groupAnnotation, 2 * g);
      this.groupItem = grow(this.groupItem, 2 * g);
      this.groupScore = grow(this.groupScore, 2 * g);
    }
    this.groupStart[g] = this.open = this.runCount;
    this.groupAnnotation[g] = annotation;
    this.groupItem[g] = item;
    this.groupScore[g] = score;
    this.sorted = true;
  }

  /** A group tagged as group g of `s` (its item and score): the score is copied between typed arrays, never boxed. */
  beginAs(s: SeqSelection, g: number): void {
    this.begin(s.groupAnnotation[g], s.groupItem[g]);
    this.groupScore[this.groupCount] = s.groupScore[g];
  }

  /** Appends a run to the open group, merged into the last when it continues it (same system, touching). */
  add(kind: number, id: number, row: number, start: number, stop: number): void {
    if (start >= stop) return;
    const k = this.runCount - 1;
    if (k >= this.open) {
      const c = cmpSys(this, k, kind, id, row);
      if (c === 0 && start >= this.start[k] && start <= this.stop[k]) {
        if (stop > this.stop[k]) this.stop[k] = stop;
        return;
      }
      if (c > 0 || (c === 0 && start < this.start[k])) this.sorted = false;
    }
    if (k + 1 === this.kind.length) {
      const self = this as unknown as Record<string, Int32Array>;
      for (const f of ["kind", "id", "row", "start", "stop"]) self[f] = grow(self[f], 0);
    }
    const i = ++this.runCount - 1;
    this.kind[i] = kind;
    this.id[i] = id;
    this.row[i] = row;
    this.start[i] = start;
    this.stop[i] = stop;
    const g = this.guard;
    if (g && (i & 1023) === 1023) {
      g.tick();
      if (this.runCount > g.maxRuns) g.over("runs", this.runCount);
    }
  }

  /** Closes the open group: sorted and merged if it came out of order; dropped when it has no runs. */
  end(): void {
    if (!this.sorted) normalize(this, this.open, this.runCount);
    if (this.runCount === this.open) return;
    this.groupStart[++this.groupCount] = this.runCount;
    if (this.guard && this.groupCount > this.guard.maxGroups) this.guard.over("groups", this.groupCount);
  }

  /** One group of one run. */
  group(kind: number, id: number, row: number, start: number, stop: number, annotation = -1, item = -1): void {
    this.begin(annotation, item);
    this.add(kind, id, row, start, stop);
    this.end();
  }

  /** The selection built: exact-size copies, the builder free to reuse. */
  build(set: SeqSet): SeqSelection {
    const [g, n] = [this.groupCount, this.runCount];
    return {
      set,
      groupCount: g,
      groupStart: this.groupStart.slice(0, g + 1),
      groupAnnotation: this.groupAnnotation.slice(0, g),
      groupItem: this.groupItem.slice(0, g),
      groupScore: this.groupScore.slice(0, g),
      runCount: n,
      kind: this.kind.slice(0, n),
      id: this.id.slice(0, n),
      row: this.row.slice(0, n),
      start: this.start.slice(0, n),
      stop: this.stop.slice(0, n),
    };
  }
}

// Sorting a group that came out of order (items over several targets, inputs): an index sort, then merged in place.
let order = new Int32Array(64);
let tmp = new Int32Array(64);
let sorting: Runs = null as never;
const byRun = (a: number, b: number) =>
  cmpSys(sorting, a, sorting.kind[b], sorting.id[b], sorting.row[b]) || sorting.start[a] - sorting.start[b];

function normalize(b: SelectionBuilder, from: number, to: number): void {
  const n = to - from;
  if (order.length < n) [order, tmp] = [new Int32Array(2 * n), new Int32Array(2 * n)];
  for (let i = 0; i < n; i++) order[i] = from + i;
  sorting = b;
  order.subarray(0, n).sort(byRun);
  for (const f of ["kind", "id", "row", "start", "stop"] as const) {
    const a = b[f];
    for (let i = 0; i < n; i++) tmp[i] = a[order[i]];
    for (let i = 0; i < n; i++) a[from + i] = tmp[i];
  }
  b.runCount = from;
  const score = b.groupScore[b.groupCount];
  b.begin(b.groupAnnotation[b.groupCount], b.groupItem[b.groupCount]);
  b.groupScore[b.groupCount] = score;
  for (let i = 0; i < n; i++)
    b.add(b.kind[from + i], b.id[from + i], b.row[from + i], b.start[from + i], b.stop[from + i]);
}

/** The first run of `u` in system (kind, id, row) ending after `start`; `u.runCount` or a later system when none. */
export function firstAfter(u: Runs & { runCount: number }, kind: number, id: number, row: number, start: number) {
  let lo = 0;
  let hi = u.runCount;
  while (lo < hi) {
    const m = (lo + hi) >> 1;
    const c = cmpSys(u, m, kind, id, row);
    if (c < 0 || (c === 0 && u.stop[m] <= start)) lo = m + 1;
    else hi = m;
  }
  return lo;
}

/** Run k of `a` intersected with `u` (one canonical group), appended to `out`'s open group. */
export function intersectRun(out: SelectionBuilder, a: Runs, k: number, u: SeqSelection): void {
  const kind = a.kind[k];
  const id = a.id[k];
  const row = a.row[k];
  const s = a.start[k];
  const t = a.stop[k];
  for (let j = firstAfter(u, kind, id, row, s); j < u.runCount && cmpSys(u, j, kind, id, row) === 0; j++) {
    if (u.start[j] >= t) break;
    out.add(kind, id, row, Math.max(s, u.start[j]), Math.min(t, u.stop[j]));
  }
}

/** Run k of `a` without `u` (one canonical group), appended to `out`'s open group. */
export function subtractRun(out: SelectionBuilder, a: Runs, k: number, u: SeqSelection): void {
  const kind = a.kind[k];
  const id = a.id[k];
  const row = a.row[k];
  const t = a.stop[k];
  let s = a.start[k];
  for (let j = firstAfter(u, kind, id, row, s); j < u.runCount && cmpSys(u, j, kind, id, row) === 0; j++) {
    if (u.start[j] >= t) break;
    out.add(kind, id, row, s, u.start[j]);
    s = Math.max(s, u.stop[j]);
  }
  out.add(kind, id, row, s, t);
}

/** Every run of `s` (or of one group) into `out`'s open group. */
export function addRuns(out: SelectionBuilder, s: SeqSelection, group = -1): void {
  // No destructuring: a call per group that is not inlined would allocate its array (QR2a's heap guards).
  const from = group < 0 ? 0 : s.groupStart[group];
  const to = group < 0 ? s.runCount : s.groupStart[group + 1];
  for (let k = from; k < to; k++) out.add(s.kind[k], s.id[k], s.row[k], s.start[k], s.stop[k]);
}

/**
 * One group per run of a one-group selection, its runs shared (seqspec:seqql#elements-and-groups: scans give one group
 * per run).
 */
export function eachRun(u: SeqSelection): SeqSelection {
  const n = u.runCount;
  const groupStart = new Int32Array(n + 1);
  for (let i = 0; i <= n; i++) groupStart[i] = i;
  // The literal `build` makes, property for property: one hidden class for every selection.
  return {
    set: u.set,
    groupCount: n,
    groupStart,
    groupAnnotation: new Int32Array(n).fill(-1),
    groupItem: new Int32Array(n).fill(-1),
    groupScore: new Float64Array(n).fill(Number.NaN),
    runCount: n,
    kind: u.kind,
    id: u.id,
    row: u.row,
    start: u.start,
    stop: u.stop,
  };
}

/** Whether runs j and k are in one system. */
export const sameSystem = (s: Runs, j: number, k: number) =>
  s.kind[j] === s.kind[k] && s.id[j] === s.id[k] && s.row[j] === s.row[k];

/** Whether two groups hold the same runs and item. */
export function sameGroup(a: SeqSelection, g: number, b: SeqSelection, h: number): boolean {
  const i = a.groupStart[g];
  const j = b.groupStart[h];
  const n = a.groupStart[g + 1] - i;
  if (n !== b.groupStart[h + 1] - j || a.groupItem[g] !== b.groupItem[h]) return false;
  if (a.groupAnnotation[g] !== b.groupAnnotation[h]) return false;
  for (let k = 0; k < n; k++) {
    if (cmpSys(a, i + k, b.kind[j + k], b.id[j + k], b.row[j + k]) !== 0) return false;
    if (a.start[i + k] !== b.start[j + k] || a.stop[i + k] !== b.stop[j + k]) return false;
  }
  return true;
}

/** A molecule run as runs on each alive unit over it, clipped to its range; false when it has none. */
export function liftMolecule(out: SelectionBuilder, set: SeqSet, m: number, s: number, t: number, units: Units) {
  const u = set.units;
  for (let k = units.start[m]; k < units.start[m + 1]; k++) {
    const x = units.units[k];
    out.add(K.unit, x, -1, Math.max(s, u.rangeStart[x]), Math.min(t, u.rangeStop[x]));
  }
  return units.start[m + 1] > units.start[m];
}
type Units = { readonly start: Int32Array; readonly units: Int32Array };
