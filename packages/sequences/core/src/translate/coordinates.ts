/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Runs in coordinate systems named by key (seqspec:libraries#coordinates): the plain-data form a host, an event or a
// command holds a target in, and its conversions to and from `Extents` and loci. `toExtents` is allocation-free once
// `out` has room (a host's hover goes through it at pointer rate); the loci conversions are edge rate and allocate.

import {
  type CoordinateSystem,
  CoordinateSystemKind as K,
  Extents,
  type Loci,
  OrderedSet,
  type SeqLoci,
  type SeqSet,
  Systems,
} from '@molstar/seq-core/model';
import type { Translator } from './translator.js';

/** Runs in coordinate systems named by key: `(index into systems, start, stop)` triples. */
export interface Coordinates {
  readonly systems: readonly string[];
  readonly runs: Int32Array;
}

/** The key prefix of a set's own systems: `seqstar:<ref>` (a molecule, a unit, a complex, columns, a row). */
export const SeqSystemPrefix = 'seqstar:';

const isCoordinates = (t: unknown): t is Coordinates =>
  typeof t === 'object' && t !== null && (t as Coordinates).runs !== undefined;

/**
 * `count` start-stop pairs (the first `2 * count` of `pairs`), sorted and merged where they touch or overlap, as an
 * `OrderedSet`: a range for one run. Edge rate: O(runs) unless several runs remain, then their positions in one array.
 */
function ofRuns(pairs: ArrayLike<number>, count: number = pairs.length / 2): OrderedSet {
  const order = Array.from({ length: count }, (_, i) => i).sort((a, b) => pairs[2 * a] - pairs[2 * b]);
  const merged: number[] = [];
  for (const i of order) {
    const [s, t] = [pairs[2 * i], pairs[2 * i + 1]];
    const n = merged.length;
    if (t <= s) continue;
    if (n > 0 && s <= merged[n - 1]) merged[n - 1] = Math.max(merged[n - 1], t);
    else merged.push(s, t);
  }
  if (merged.length === 2) return OrderedSet.ofRange(merged[0], merged[1]);
  let total = 0;
  for (let i = 0; i < merged.length; i += 2) total += merged[i + 1] - merged[i];
  const out = new Int32Array(total);
  let at = 0;
  for (let i = 0; i < merged.length; i += 2) for (let x = merged[i]; x < merged[i + 1]; x++) out[at++] = x;
  return OrderedSet.ofSortedArray(out);
}

/**
 * `Coordinates` appended to `out` as entries in the set's own systems, through `t` (the second step of every
 * conversion): a key `t` does not know (or that is unregistered) gives nothing; a host's external system's runs go
 * through each map registered with it into the system at its other end; a complex or an item system gives nothing (a
 * complex names its members). `scratch`: a one-entry `Extents` the call reuses. The entries added; allocation-free once
 * `out` has room.
 */
function toExtents(t: Translator, c: Coordinates, out: Extents, scratch: Extents): number {
  const before = out.count;
  const { runs, systems } = c;
  for (let r = 0; r + 2 < runs.length; r += 3) {
    const key = systems[runs[r]];
    const s = t.system(key);
    if (!s || !t.descriptor(key)) continue;
    if (s.kind === K.external) external(t, s, runs[r + 1], runs[r + 2], out, scratch);
    else if (s.kind !== K.complex && s.kind !== K.item) Extents.add(out, s, runs[r + 1], runs[r + 2]);
  }
  return out.count - before;
}

/** An external run through each map registered with its system, into the Seq* system at its other end. */
function external(t: Translator, s: CoordinateSystem, start: number, stop: number, out: Extents, one: Extents) {
  const near = t.neighbors(s);
  for (let k = 0; k < near.length; k++) {
    const n = near[k];
    if (n.kind === K.external || n.kind === K.complex) continue;
    one.count = 0;
    Extents.add(one, s, start, stop);
    t.translate(one, s, n, out);
  }
}

/**
 * Loci as `Coordinates` in their set's own systems (edge rate, allocates): each element's runs under its system's key,
 * an item as its rows' ranges; a unit of a complex also lists the complex, whose descriptor's `members` name its units,
 * so a host builds one structure loci from the runs.
 */
function ofLoci(loci: Loci): Coordinates {
  const ext = Extents.create(8);
  Extents.addLoci(ext, loci);
  return ofExtents(loci.set, ext);
}

/** Entries in `set`'s own systems as `Coordinates` (edge rate, allocates); a kind no key names is left out. */
function ofExtents(set: SeqSet, ext: Extents, keyOf?: (s: CoordinateSystem) => string | null): Coordinates {
  const systems: string[] = [];
  const runs: number[] = [];
  const index = (key: string) => {
    const k = systems.indexOf(key);
    return k >= 0 ? k : systems.push(key) - 1;
  };
  const sys = { kind: 0, id: 0, row: -1 };
  for (let i = 0; i < ext.count; i++) {
    sys.kind = ext.kind[i];
    sys.id = ext.id[i];
    sys.row = ext.row[i];
    const ref = sys.kind === K.external ? null : Systems.refOf(set, sys.kind, sys.id, sys.row);
    const key = ref !== null ? SeqSystemPrefix + ref : (keyOf?.(sys) ?? null);
    if (key !== null) runs.push(index(key), ext.start[i], ext.stop[i]);
    const cx = ext.kind[i] === K.unit ? set.complexes[set.units.complex[ext.id[i]]] : null;
    if (cx) index(SeqSystemPrefix + set.strings.get(cx.key));
  }
  return { systems, runs: Int32Array.from(runs) };
}

/**
 * The unit and molecule entries of `ext` as one seq loci of `set`, a molecule entry in every alive unit of its
 * molecule (every copy); null: none. Edge rate: per unit its runs, never a Set of positions, so a chromosome-long run
 * converts in O(runs).
 */
function toSeqLoci(set: SeqSet, ext: Extents): SeqLoci | null {
  const byUnit = new Map<number, number[]>();
  const put = (u: number, s: number, t: number) => {
    const list = byUnit.get(u) ?? byUnit.set(u, []).get(u);
    list?.push(s, t);
  };
  for (let i = 0; i < ext.count; i++) {
    if (ext.kind[i] === K.unit) put(ext.id[i], ext.start[i], ext.stop[i]);
    if (ext.kind[i] !== K.molecule) continue;
    for (let u = 0; u < set.units.count; u++)
      if (set.units.alive[u] && set.units.molecule[u] === ext.id[i]) put(u, ext.start[i], ext.stop[i]);
  }
  const elements = [...byUnit].sort((a, b) => a[0] - b[0]).map(([unit, runs]) => ({ unit, indices: ofRuns(runs) }));
  return elements.length ? ({ kind: 'seq-loci', set, elements } as unknown as SeqLoci) : null;
}

/**
 * Entries in `set`'s systems as loci (edge rate, allocates): its units and molecules as one seq loci, rows as one row
 * loci and columns as one column loci per alignment. Other kinds are left out.
 */
function toLoci(set: SeqSet, ext: Extents): Loci[] {
  const seq = toSeqLoci(set, ext);
  const out: Loci[] = seq ? [seq] : [];
  const runs = new Map<number, Map<number, number[]>>(); // alignment → row (-1: its columns) → runs
  for (let i = 0; i < ext.count; i++) {
    if (ext.kind[i] !== K.row && ext.kind[i] !== K.columns) continue;
    const byRow = runs.get(ext.id[i]) ?? runs.set(ext.id[i], new Map()).get(ext.id[i]);
    const row = ext.kind[i] === K.row ? ext.row[i] : -1;
    byRow?.set(row, [...(byRow.get(row) ?? []), ext.start[i], ext.stop[i]]);
  }
  for (const [alignment, byRow] of runs) {
    const rows = [...byRow].filter(([r]) => r >= 0).sort((a, b) => a[0] - b[0]);
    const elements = rows.map(([row, r]) => ({ row, indices: ofRuns(r) }));
    if (elements.length) out.push({ kind: 'row-loci', set, alignment, elements } as unknown as Loci);
    const cols = byRow.get(-1);
    if (cols) out.push({ kind: 'column-loci', set, alignment, columns: ofRuns(cols) } as unknown as Loci);
  }
  return out;
}

/** `Coordinates`' conversions (seqspec:libraries#coordinates). */
export const Coordinates = {
  is: isCoordinates,
  ofRuns,
  toExtents,
  ofLoci,
  ofExtents,
  toSeqLoci,
  toLoci,
};
