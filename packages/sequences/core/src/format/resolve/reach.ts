/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// What a section reaches (seqspec:model#coordinate-systems, seqspec:format#placement): an annotation's ends (its side-a
// slots, and a pair's side-b targets), paths from an end's coordinate system to a block, and the positions an end
// covers there. Shared by `positions`, focus, placement and links. Resolve-time code: it allocates per call, never per
// position, and the hover path never runs it.

import {
  type Annotation,
  Annotations,
  type ComposedPath,
  type CoordinateSystem,
  CoordinateSystemKind,
  Mappings,
  type PathStep,
  Paths,
  type SeqSet,
  SortedArray,
  type Targets,
} from '@molstar/seq-core/model';
import { numberColumn } from '../palettes/fields.js';
import type { BlockRow } from './section.js';

/** A side-a slot, or a distinct side-b target of a pairs annotation (slot = slot count + its first-row order). */
export interface End {
  readonly system: CoordinateSystem;
  readonly slot: number;
  readonly side: 0 | 1;
  readonly rows: number[];
}

export const unitSystem = (u: number): CoordinateSystem => ({ kind: CoordinateSystemKind.unit, id: u, row: -1 });
export const columnsSystem = (a: number): CoordinateSystem => ({ kind: CoordinateSystemKind.columns, id: a, row: -1 });
export const blockSystem = (b: { alignment: number; unit: number }) =>
  b.alignment >= 0 ? columnsSystem(b.alignment) : unitSystem(b.unit);

function systemAt(t: Targets, i: number): CoordinateSystem | null {
  if (t.unit[i] >= 0) return unitSystem(t.unit[i]);
  if (t.molecule[i] >= 0) return { kind: CoordinateSystemKind.molecule, id: t.molecule[i], row: -1 };
  if (t.alignment[i] < 0) return null;
  const row = t.row[i];
  return row >= 0 ? { kind: CoordinateSystemKind.row, id: t.alignment[i], row } : columnsSystem(t.alignment[i]);
}

const ends = new WeakMap<Annotation, End[]>();
export function endsOf(ann: Annotation): End[] {
  let out = ends.get(ann);
  if (out) return out;
  out = [];
  const bySlot: (End | null)[] = [];
  const n = ann.slots.unit.length;
  for (let s = 0; s < n; s++) {
    const system = systemAt(ann.slots, s);
    bySlot.push(system ? { system, slot: s, side: 0, rows: [] } : null);
  }
  for (let i = 0; i < ann.count; i++) bySlot[ann.slot[i]]?.rows.push(i);
  for (const e of bySlot) if (e) out.push(e);
  const b = new Map<string, End>();
  for (let i = 0; ann.pairs && i < ann.count; i++) {
    const system = systemAt(ann.pairs.bTarget, i);
    if (!system) continue;
    const k = `${system.kind}:${system.id}:${system.row}`;
    const e = b.get(k) ?? { system, slot: n + b.size, side: 1, rows: [] };
    b.set(k, e);
    e.rows.push(i);
  }
  for (const e of b.values()) out.push(e);
  ends.set(ann, out);
  return out;
}

/** The last path found (`pathTo`): its steps. */
export const steps: PathStep[] = [];

let searched = 0;
/** The path searches `pathTo` ran (V7's counter: a rerun of `composeView` runs none). */
export const pathsSearched = () => searched;

/**
 * The composed path from one coordinate system to another (null: the same; undefined: none). Copies never reach each
 * other: the model's search refuses a path up from a unit and down to another unit of its molecule.
 */
export function pathTo(
  set: SeqSet,
  from: CoordinateSystem,
  to: CoordinateSystem,
): ComposedPath | null | undefined | 'ambiguous' {
  searched++;
  const n = Mappings.path(set, from, to, steps);
  if (n === -1 || n === 0) return n === 0 ? null : undefined;
  return n === -2 ? 'ambiguous' : Mappings.compose(set, from, steps);
}

/** The row whose edge the last path crossed (its lane), or -1. */
const rowOf = () => {
  for (let k = steps.length - 1; k >= 0; k--) if (steps[k].mapping < 0 && steps[k].unit < 0) return steps[k].row;
  return -1;
};

/**
 * Calls back per (member, path, row) by which a coordinate system reaches a block: a molecule every unit over it alike
 * (member -1), a unit itself; anything else along a path. A merged block's member is the unit landed on. False:
 * ambiguous.
 */
export function reach(
  set: SeqSet,
  from: CoordinateSystem,
  block: BlockRow,
  cb: (member: number, path: ComposedPath | null, row: number) => void,
): boolean {
  if (block.alignment >= 0) {
    const p = pathTo(set, from, columnsSystem(block.alignment));
    if (p !== undefined && p !== 'ambiguous') cb(-1, p, p ? rowOf() : -1);
    return p !== 'ambiguous';
  }
  if (from.kind === CoordinateSystemKind.molecule && from.id === set.units.molecule[block.unit]) {
    cb(-1, null, -1);
    return true;
  }
  const merged = block.members !== null;
  // A unit onto its own copies: only itself, unless a mapping ends on a unit (seqspec:model#copies); no search per copy.
  if (
    from.kind === CoordinateSystemKind.unit &&
    set.units.molecule[from.id] === set.units.molecule[block.unit] &&
    !Mappings.unitsMapped(set)
  ) {
    if (block.members ? SortedArray.has(block.members, from.id) : from.id === block.unit)
      cb(merged ? from.id : -1, null, -1);
    return true;
  }
  let ok = true;
  for (const m of block.members ?? [block.unit]) {
    if (from.kind === CoordinateSystemKind.unit && from.id === m) {
      cb(merged ? m : -1, null, -1);
      continue;
    }
    const p = pathTo(set, from, unitSystem(m));
    ok &&= p !== 'ambiguous';
    if (!p || p === 'ambiguous') continue;
    const last = steps[steps.length - 1];
    const alike = last.unit >= 0 && last.reversed;
    cb(merged && !alike ? m : -1, p, rowOf());
    if (alike) break;
  }
  return ok;
}

/**
 * Marks the positions an end covers (seqspec:format#sections): ranges their rows, pairs their endpoints, values and
 * matrix the positions with a value; through `path` into the block's coordinate system (`mask`, by its indices).
 * `keep` filters rows. Nothing is allocated per row or position (resolve.heap.ts).
 */
export function cover(
  ann: Annotation,
  end: End,
  path: ComposedPath | null,
  mask: Uint8Array,
  keep?: (row: number) => boolean,
) {
  const dense = end.side === 0 ? ann.dense?.[end.slot] : undefined;
  if (dense) for (let i = 0; i < dense.length; i++) if (!Number.isNaN(dense[i])) mark(path, i, mask);
  // A row with no value covers nothing only in `values` geometry; a ranges row's `value` is optional.
  const valued = ann.geometry === 'values' ? ann.sparse : null;
  const rows = end.rows;
  for (let k = 0; k < rows.length; k++) {
    const r = rows[k];
    if (keep && !keep(r)) continue;
    if (end.side === 1 && ann.pairs) {
      const b = ann.pairs.bStart[r];
      if (b < 0) span(path, 0, path ? path.length : mask.length, mask);
      else span(path, b, ann.pairs.bStop[r], mask);
    } else if (!(valued && Number.isNaN(valued[r]))) span(path, ann.start[r], ann.stop[r], mask);
  }
}
function span(path: ComposedPath | null, s: number, t: number, mask: Uint8Array) {
  for (let i = s; i < Math.max(t, s + 1); i++) mark(path, i, mask);
}
/** Position i's image; through a codon step from the protein, its codon's three (seqspec:model#codon-maps). */
function mark(path: ComposedPath | null, i: number, mask: Uint8Array) {
  if (path?.heldSpans) for (let k = 0; k < 3; k++) at(mask, Paths.spanAt(path, i, k));
  else at(mask, path ? Paths.at(path, i) : i);
}
function at(mask: Uint8Array, j: number) {
  if (j >= 0 && j < mask.length) mask[j] = 1;
}

/**
 * A filter's test of a row's field (seqspec:format#layer-types: values compared as strings), allocating nothing per
 * row: a number column is compared as numbers (a double read through `rowValue` would be boxed), anything else by its
 * string.
 */
export function matcher(values: readonly unknown[]) {
  const strings = new Set(values.map(String));
  const numbers = Float64Array.from(values.filter((v) => String(Number(v)) === String(v)).map(Number));
  return (a: Annotation, set: SeqSet, row: number, field: string) => {
    const column = numberColumn(a, set, field);
    if (column) {
      for (let i = 0; i < numbers.length; i++) if (numbers[i] === column[row]) return true;
      return false;
    }
    const x = Annotations.rowValue(a, set, row, field);
    return x !== undefined && strings.has(String(x));
  };
}
