/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The sequence selection (seqspec:seqql#the-selection): groups of runs in coordinate systems, in typed arrays, each
// group optionally an annotation item; immutable once built (kernels.ts builds it). Its operations are edges: a host's
// input in, loci, extents and a summary out.

import { type AnnotationInput, Annotations } from "../annotation.ts";
import type { AnnotationId, ItemIndex, UnitId } from "../common.ts";
import { CoordinateSystemKind as K } from "../coordinates.ts";
import { Extents } from "../extents.ts";
import type { ColumnLoci, ItemLoci, Loci, RowLoci, SeqLoci } from "../loci.ts";
import { OrderedSet } from "../ordered-set.ts";
import type { SeqSet } from "../set.ts";
import { addRuns, liftMolecule, SelectionBuilder, sameSystem } from "./kernels.ts";
import { domainOf, labelAt, nameOf, numberingOf, unitsOf } from "./systems.ts";

export interface SeqSelection {
  readonly set: SeqSet;
  readonly groupCount: number;
  /** groupCount + 1 offsets into the runs. */
  readonly groupStart: Int32Array;
  /** AnnotationId of the group's item, -1: none. */
  readonly groupAnnotation: Int32Array;
  /** ItemIndex, -1: none. */
  readonly groupItem: Int32Array;
  /** A hit's score, NaN: none. */
  readonly groupScore: Float64Array;
  readonly runCount: number;
  /** CoordinateSystemKind: molecule, unit, row, columns. */
  readonly kind: Uint8Array;
  readonly id: Int32Array;
  readonly row: Int32Array;
  /** Half-open, start < stop. */
  readonly start: Int32Array;
  readonly stop: Int32Array;
}

const domain = new Int32Array(2);
const scratch = Extents.create(64);

/** One group per entry of `groups`, each an `Extents` slice: entries clipped to their system's domain. */
function fromEntries(set: SeqSet, e: Extents, cuts: readonly number[]): SeqSelection {
  const b = new SelectionBuilder();
  for (let g = 0; g + 1 < cuts.length; g++) {
    b.begin();
    for (let i = cuts[g]; i < cuts[g + 1]; i++) {
      if (!domainOf(set, e.kind[i], e.id[i], e.row[i], domain)) continue;
      b.add(e.kind[i], e.id[i], e.row[i], Math.max(e.start[i], domain[0]), Math.min(e.stop[i], domain[1]));
    }
    b.end();
  }
  return b.build(set);
}

/** One canonical group of every run of `s` (a group of an item keeps nothing of it). */
function union(s: SeqSelection): SeqSelection {
  if (s.groupCount === 1 && s.groupAnnotation[0] < 0) return s;
  const b = new SelectionBuilder();
  b.begin();
  addRuns(b, s);
  b.end();
  return b.build(s.set);
}

function positionCount(s: SeqSelection, group?: number): number {
  const u = group === undefined ? union(s) : s;
  let n = 0;
  const [from, to] = group === undefined ? [0, u.runCount] : [s.groupStart[group], s.groupStart[group + 1]];
  for (let k = from; k < to; k++) n += u.stop[k] - u.start[k];
  return n;
}

/** Runs [k, j) of one system as an ordered set: an Interval for one run, else its positions. */
function indicesOf(s: SeqSelection, k: number, j: number): OrderedSet {
  if (j - k === 1) return OrderedSet.ofRange(s.start[k], s.stop[k]);
  let n = 0;
  for (let i = k; i < j; i++) n += s.stop[i] - s.start[i];
  const xs = new Int32Array(n);
  for (let i = k, at = 0; i < j; i++) for (let x = s.start[i]; x < s.stop[i]; x++) xs[at++] = x;
  return OrderedSet.ofSortedArray(xs);
}

export interface ToLociOptions {
  /** One group's loci; default every group's. */
  readonly group?: number;
  /** "positions": tagged groups as their runs, not as items. */
  readonly items?: "items" | "positions";
  /** Receives the positions of molecule runs dropped for want of an alive unit over the molecule. */
  readonly report?: { dropped: number };
}

/**
 * The model's loci, at most one per loci kind and alignment: unit runs a `SeqLoci`, molecule runs that range
 * in every alive unit over the molecule (dropped when none), row runs a `RowLoci`, column runs a `ColumnLoci`, tagged
 * groups an `ItemLoci` (their runs instead with `items: "positions"`). O(positions) for multi-run systems.
 */
function toLoci(s: SeqSelection, opts: ToLociOptions = {}): readonly Loci[] {
  const { set } = s;
  const units = unitsOf(set);
  const items = new Map<number, number[]>();
  const b = new SelectionBuilder();
  let dropped = 0;
  b.begin();
  const [g0, g1] = opts.group === undefined ? [0, s.groupCount] : [opts.group, opts.group + 1];
  for (let g = g0; g < g1; g++) {
    const a = s.groupAnnotation[g];
    if (a >= 0 && opts.items !== "positions") {
      const list = items.get(a) ?? [];
      list.push(s.groupItem[g]);
      items.set(a, list);
      continue;
    }
    for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
      if (s.kind[k] !== K.molecule) b.add(s.kind[k], s.id[k], s.row[k], s.start[k], s.stop[k]);
      else if (!liftMolecule(b, set, s.id[k], s.start[k], s.stop[k], units)) dropped += s.stop[k] - s.start[k];
    }
  }
  if (opts.report) opts.report.dropped = dropped;
  b.end();
  const u = b.build(set);
  const out: Loci[] = [];
  const seq: { unit: UnitId; indices: OrderedSet }[] = [];
  let rows: { row: number; indices: OrderedSet }[] = [];
  for (let k = 0, j = 0; k < u.runCount; k = j) {
    while (j < u.runCount && sameSystem(u, j, k)) j++;
    const [kind, id, indices] = [u.kind[k], u.id[k], indicesOf(u, k, j)];
    if (kind === K.unit) seq.push({ unit: id as UnitId, indices });
    else if (kind === K.columns) out.push({ kind: "column-loci", set, alignment: id, columns: indices } as ColumnLoci);
    else if (kind === K.row) {
      const last = out.at(-1);
      if (last?.kind !== "row-loci" || last.alignment !== id) {
        rows = [];
        out.push({ kind: "row-loci", set, alignment: id, elements: rows } as unknown as RowLoci);
      }
      rows.push({ row: u.row[k], indices });
    }
  }
  if (seq.length) out.unshift({ kind: "seq-loci", set, elements: seq } as SeqLoci);
  const elements = [...items.keys()]
    .sort((x, y) => x - y)
    .map((a) => {
      const sorted = Int32Array.from(new Set(items.get(a))).sort();
      return { annotation: a as AnnotationId, items: OrderedSet.ofSortedArray<ItemIndex>(sorted) };
    });
  if (elements.length) out.push({ kind: "item-loci", set, elements } as ItemLoci);
  return out;
}

/** Appends the runs (of one group) to `out`; the entries added. Allocates nothing once `out` has room. */
function toExtents(s: SeqSelection, out: Extents, group?: number): number {
  const [from, to] = group === undefined ? [0, s.runCount] : [s.groupStart[group], s.groupStart[group + 1]];
  const sys = { kind: 0, id: 0, row: 0 };
  for (let k = from; k < to; k++) {
    sys.kind = s.kind[k];
    sys.id = s.id[k];
    sys.row = s.row[k];
    Extents.add(out, sys, s.start[k], s.stop[k]);
  }
  return to - from;
}

const count = (n: number, what: string) => `${n.toLocaleString("en-US")} ${what}${n === 1 ? "" : "s"}`;

/**
 * The line a person or the agent reads: "3 groups, 412 positions: 1brs/A 10-20, 45-60; msa columns 101-180
 * (and 1 more)", in labels of each system's default numbering; `indices` adds each run's slice (`10-20 @9:20`).
 */
function summary(s: SeqSelection, opts: { maxParts?: number; indices?: boolean } = {}): string {
  const u = union(s);
  if (u.runCount === 0) return "selects nothing";
  const { set } = s;
  const max = opts.maxParts ?? 8;
  const parts: string[] = [];
  let n = null as ReturnType<typeof numberingOf>;
  const shown = Math.min(max, u.runCount);
  for (let k = 0; k < shown; k++) {
    const fresh = k === 0 || !sameSystem(u, k - 1, k);
    if (fresh) n = numberingOf(set, u.kind[k], u.id[k], u.row[k], null);
    const [a, z] = [labelAt(set, n ?? null, u.start[k]), labelAt(set, n ?? null, u.stop[k] - 1)];
    // An endpoint without a label (a gap in the numbering): the run by its slice only.
    const slice = ` @${u.start[k]}:${u.stop[k]}`;
    const text = !a || !z ? slice.trim() : `${a === z ? a : `${a}-${z}`}${opts.indices ? slice : ""}`;
    if (fresh) parts.push(`${nameOf(set, u.kind[k], u.id[k], u.row[k])} ${text}`);
    else parts[parts.length - 1] += `, ${text}`;
  }
  const more = u.runCount > shown ? ` (and ${u.runCount - shown} more)` : "";
  return `${count(s.groupCount, "group")}, ${count(positionCount(u), "position")}: ${parts.join("; ")}${more}`;
}

/**
 * A run's system as a `ranges` target (seqspec:format#annotation-tables): a unit's or molecule's key,
 * `<alignment>/<row>`, an alignment's.
 */
function targetOf(set: SeqSet, kind: number, id: number, row: number): string {
  if (kind === K.unit) return set.strings.get(set.units.key[id]);
  if (kind === K.molecule) return set.strings.get(set.molecules[id]?.key ?? -1);
  const a = set.strings.get(set.alignments[id]?.key ?? -1);
  return kind === K.row ? `${a}/${row}` : a;
}

export interface RangesInputOptions {
  /** The annotation's label. */
  readonly label?: string;
  /** Per group, its item's id; default the group's item's id, else `q<g>`. */
  readonly ids?: readonly string[];
  /** Per group, its rows' label; default the group's item's label. */
  readonly labels?: readonly (string | null)[];
}
/**
 * A location set (LS1, LS2): a `ranges` annotation `key`, one item per group, its runs as rows on their systems' refs,
 * positions as indices; `label` and `score` columns (the group's score, NaN: none). Edge rate.
 */
function toRangesInput(s: SeqSelection, key: string, opts: RangesInputOptions = {}): AnnotationInput {
  const set = s.set;
  const c = { target: [] as string[], start: [] as number[], stop: [] as number[], id: [] as string[] };
  const more = { label: [] as (string | null)[], score: [] as number[] };
  for (let g = 0; g < s.groupCount; g++) {
    const a = set.annotations[s.groupAnnotation[g]];
    const item = a ? s.groupItem[g] : -1;
    const own = a && item >= 0 && a.itemKey[item] >= 0 ? set.strings.get(a.itemKey[item]) : null;
    const rows = a && item >= 0 ? Annotations.itemRows(a) : null;
    const first = a && rows && rows.start[item] < rows.start[item + 1] ? rows.rows[rows.start[item]] : -1;
    const itemLabel = a && first >= 0 ? Annotations.rowValue(a, set, first, "label") : undefined;
    const label = opts.labels?.[g] ?? (typeof itemLabel === "string" ? itemLabel : null);
    for (let k = s.groupStart[g]; k < s.groupStart[g + 1]; k++) {
      c.target.push(targetOf(set, s.kind[k], s.id[k], s.row[k]));
      c.start.push(s.start[k]);
      c.stop.push(s.stop[k]);
      c.id.push(opts.ids?.[g] ?? own ?? `q${g}`);
      more.label.push(label);
      more.score.push(s.groupScore[g]);
    }
  }
  return { key, label: opts.label, geometry: "ranges", columns: { ...c, ...more } };
}

export const SeqSelection = {
  empty: (set: SeqSet): SeqSelection => new SelectionBuilder().build(set),
  /** A host's input: one group per loci, through `Extents.addLoci` (boundaries are not elements, and drop out). */
  fromLoci(set: SeqSet, loci: readonly Loci[]): SeqSelection {
    scratch.count = 0;
    const cuts = [0];
    for (const l of loci) cuts.push(cuts[cuts.length - 1] + Extents.addLoci(scratch, l));
    return fromEntries(set, scratch, cuts);
  },
  /** One group. */
  fromExtents: (set: SeqSet, e: Extents): SeqSelection => fromEntries(set, e, [0, e.count]),
  isEmpty: (s: SeqSelection) => s.groupCount === 0,
  /** Of a group; else of the union (positions in several groups counted once). */
  positionCount,
  /** Allocation-free: `visit(k, c)` per run index of the group. */
  forEachRun<C>(s: SeqSelection, group: number, visit: (k: number, c: C) => void, c: C): void {
    for (let k = s.groupStart[group]; k < s.groupStart[group + 1]; k++) visit(k, c);
  },
  union,
  toExtents,
  toLoci,
  toRangesInput,
  summary,
};
