/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Runs in coordinate systems, reused (seqspec:model#extents; seqspec:libraries#translate): entries of (system,
// start, stop), half-open, which may repeat. Beside `CoordinateMap` since VW5 (TR1 had it in translate): the projection
// core (the format's `resolve/project.ts`) and `translate` read this one declaration.

import { Annotations } from './annotation.js';
import { type CoordinateSystem, CoordinateSystemKind as K } from './coordinates.js';
import type { ItemLoci, Loci } from './loci.js';
import { OrderedSet } from './ordered-set.js';
import type { SeqSet } from './set.js';

export interface Extents {
  count: number;
  kind: Uint8Array;
  id: Int32Array;
  row: Int32Array;
  start: Int32Array;
  stop: Int32Array;
}

function create(capacity = 16): Extents {
  const n = Math.max(1, capacity);
  const [id, row, start, stop] = [0, 1, 2, 3].map(() => new Int32Array(n));
  return { count: 0, kind: new Uint8Array(n), id, row, start, stop };
}

/** Appends one run; the arrays double when full (the first uses only, for a reused `Extents`). */
function add(e: Extents, s: CoordinateSystem, start: number, stop: number): void {
  if (e.count === e.kind.length) grow(e);
  const i = e.count++;
  e.kind[i] = s.kind;
  e.id[i] = s.id;
  e.row[i] = s.row;
  e.start[i] = start;
  e.stop[i] = stop;
}

function grow(e: Extents) {
  const n = 2 * e.kind.length;
  const kind = new Uint8Array(n);
  kind.set(e.kind);
  e.kind = kind;
  for (const f of ['id', 'row', 'start', 'stop'] as const) {
    const a = new Int32Array(n);
    a.set(e[f]);
    e[f] = a;
  }
}

/** The one scratch system `addLoci` passes to `add` (which copies it): no object per entry. */
const sys: CoordinateSystem = { kind: 0, id: -1, row: -1 };
const at = (kind: number, id: number, row = -1): CoordinateSystem => {
  sys.kind = kind;
  sys.id = id;
  sys.row = row;
  return sys;
};

/** An ordered set as runs of consecutive indices in `s` (an Interval one run; a SortedArray walked once). */
function addRuns(out: Extents, s: CoordinateSystem, set: OrderedSet): void {
  if (typeof set === 'number') {
    if (OrderedSet.size(set) > 0) add(out, s, OrderedSet.min(set), OrderedSet.max(set) + 1);
    return;
  }
  for (let i = 0; i < set.length; ) {
    let j = i;
    while (j + 1 < set.length && set[j + 1] === set[j] + 1) j++;
    add(out, s, set[i], set[j] + 1);
    i = j + 1;
  }
}

/** An aligned row's columns: its residues through its map in the row's system, a gap the column in the columns'. */
function addCells(out: Extents, set: SeqSet, alignment: number, row: number, columns: OrderedSet): void {
  const rows = set.alignments[alignment]?.rows;
  let mode = 0; // 1: residues of the row, 2: columns
  let a = 0;
  let b = 0;
  for (let i = 0, n = OrderedSet.size(columns); i < n; i++) {
    const c = OrderedSet.getAt(columns, i);
    const p = rows ? rows.colToPos(row as never, c as never) : -1;
    const m = p >= 0 ? 1 : 2;
    const x = p >= 0 ? p : c;
    if (m === mode && x === b) b++;
    else {
      if (mode > 0) add(out, mode === 1 ? at(K.row, alignment, row) : at(K.columns, alignment), a, b);
      mode = m;
      a = x;
      b = x + 1;
    }
  }
  if (mode > 0) add(out, mode === 1 ? at(K.row, alignment, row) : at(K.columns, alignment), a, b);
}

/** An item row's target (a slot of `t`, or a pair's side-b row) as a system in the scratch; false: none. */
function target(t: { unit: Int32Array; molecule: Int32Array; alignment: Int32Array; row: Int32Array }, i: number) {
  if (t.unit[i] >= 0) at(K.unit, t.unit[i]);
  else if (t.molecule[i] >= 0) at(K.molecule, t.molecule[i]);
  else if (t.alignment[i] >= 0) at(t.row[i] >= 0 ? K.row : K.columns, t.alignment[i], t.row[i]);
  else return false;
  return true;
}

/**
 * Appends `loci` to `out` (not cleared) as entries in each element's own coordinate system: runs of consecutive
 * indices (an `Interval` one run; a `SortedArray` walked once). Returns the entries added. Allocates nothing once `out`
 * has room. A boundary is `(b, b)`; an aligned cell is its residue in the row's system, a gap its column; an item is
 * its rows' ranges in their targets' systems (a pair's side b too).
 */
function addLoci(out: Extents, loci: Loci): number {
  const before = out.count;
  const set = loci.set;
  if (loci.kind === 'seq-loci') {
    for (const e of loci.elements) addRuns(out, at(K.unit, e.unit), e.indices);
  } else if (loci.kind === 'boundary-loci') {
    for (const e of loci.elements)
      for (let i = 0, n = OrderedSet.size(e.indices); i < n; i++) {
        const x = OrderedSet.getAt(e.indices, i);
        add(out, at(K.unit, e.unit), x, x);
      }
  } else if (loci.kind === 'row-loci') {
    for (const e of loci.elements) addRuns(out, at(K.row, loci.alignment, e.row), e.indices);
  } else if (loci.kind === 'column-loci') {
    addRuns(out, at(K.columns, loci.alignment), loci.columns);
  } else if (loci.kind === 'aligned-loci') {
    const el = loci.elements;
    if (!('rows' in el)) for (const e of el) addCells(out, set, loci.alignment, e.row, e.columns);
    else
      for (let i = 0, n = OrderedSet.size(el.rows); i < n; i++)
        addCells(out, set, loci.alignment, OrderedSet.getAt(el.rows, i), el.columns);
  } else {
    for (const e of loci.elements) addItems(out, set, e);
  }
  return out.count - before;
}

/**
 * An item loci element's items, each as `addItem`. A function of its own: inline in `addLoci` beside the other kinds,
 * V8 13.6 stopped inlining `size` and `getAt` there and boxed the packed `Interval` into each call (2026-10-13).
 */
function addItems(out: Extents, set: SeqSet, e: ItemLoci['elements'][number]): void {
  const items = e.items;
  for (let i = 0, n = OrderedSet.size(items); i < n; i++) addItem(out, set, e.annotation, OrderedSet.getAt(items, i));
}

/** One item's rows (`itemRows`) appended to `out` as `addLoci` adds item loci: a hovered item, allocation-free. */
function addItem(out: Extents, set: SeqSet, annotation: number, item: number): void {
  const ann = set.annotations[annotation];
  if (!ann || ann.geometry === 'rows' || item < 0 || item >= ann.itemCount) return;
  const byItem = Annotations.itemRows(ann);
  for (let k = byItem.start[item]; k < byItem.start[item + 1]; k++) {
    const r = byItem.rows[k];
    if (target(ann.slots, ann.slot[r])) add(out, sys, ann.start[r], ann.stop[r]);
    if (ann.pairs && target(ann.pairs.bTarget, r)) {
      const whole = ann.pairs.bStart[r] < 0;
      add(out, sys, whole ? 0 : ann.pairs.bStart[r], whole ? 0x7fffffff : ann.pairs.bStop[r]);
    }
  }
}

/** A reused `Extents` starts again with `count = 0`. */
export const Extents = { create, add, addLoci, addItem };
