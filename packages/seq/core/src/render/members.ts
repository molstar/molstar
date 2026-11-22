/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// A categorical axis's members in order (_spec/README.md "Render", The grid layer: categorical axes; EC-35): the
// permutation from a display offset to a member (`order`) and its inverse (`rank`), both `Int32Array`, the members'
// labels, and a version that moves with every order made (a grid's pyramid in display order is keyed by it). The order
// is the rows layer's values: `data`; `tree` (a `tree` node of the alignment, its leaves matched to the rows by label,
// rows the tree lacks after them in data order); `{ field, direction }` over the rows table (`label`, an inline row
// field): numbers by value, text by code points, missing last, ties in data order. Made once a node and order, kept.

import { CategoryColumn, exactAddress, isColumn, type Node, type Value } from '@molstar/core/repr/tree';
import { readNewick } from '@molstar/seqx-core/data/alignment-readers';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import type { MemberAxis, MemberOrder } from '@molstar/seqx-core/view';

let Versions = 0;

/** A categorical axis's members in display order. */
export class Members {
  /** Each member's display offset. */
  readonly rank: Int32Array;
  /** The axis's pieces: one a member, in display order (`[m, m + 1, …]`). */
  readonly pieces: readonly number[];
  readonly version = ++Versions;
  constructor(
    /** The member at each display offset. */
    readonly order: Int32Array,
    /** Each member's label (data order). */
    readonly labels: readonly string[],
  ) {
    const n = order.length;
    this.rank = new Int32Array(n);
    const pieces = new Array<number>(2 * n);
    for (let d = 0; d < n; d++) {
      this.rank[order[d]] = d;
      pieces[2 * d] = order[d];
      pieces[2 * d + 1] = order[d] + 1;
    }
    this.pieces = pieces;
  }
}

/** A node of the spec by its exact address. */
function nodeAt(spec: Node, address: string): Node | undefined {
  let n: Node | undefined = spec;
  for (const id of address.split('/').filter(Boolean)) n = n?.children?.find((c) => c.id === id);
  return n;
}

/** A column as a list (a category column's values, a typed array's numbers, a list). */
function listOf(v: Value | undefined): readonly unknown[] | undefined {
  if (v instanceof CategoryColumn) return Array.from(v.codes, (k) => (k < 0 ? null : v.dict[k]));
  if (isColumn(v)) return Array.from(v as ArrayLike<number>);
  return Array.isArray(v) ? v : undefined;
}

/** The members' labels: an inline alignment's row labels, an inline matrix's states; else `row 1`, `state 1`, …. */
function labelsOf(owner: Node | undefined, axis: MemberAxis, n: number): string[] {
  const rows = owner?.params?.rows;
  const given =
    axis.kind === 'states'
      ? listOf(owner?.params?.states)
      : rows && typeof rows === 'object' && !Array.isArray(rows) && !isColumn(rows)
        ? listOf((rows as Record<string, Value>).label)
        : undefined;
  const word = axis.kind === 'states' ? 'state' : 'row';
  const out: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const v = given?.[i];
    out[i] = v === null || v === undefined ? `${word} ${i + 1}` : String(v);
  }
  return out;
}

/** A field's values a member: `label`, else an inline row field (missing: undefined). */
function fieldOf(owner: Node | undefined, field: string, labels: readonly string[]): readonly unknown[] {
  if (field === 'label') return labels;
  const rows = owner?.params?.rows;
  const col =
    rows && typeof rows === 'object' && !Array.isArray(rows) ? (rows as Record<string, Value>)[field] : undefined;
  return listOf(col) ?? [];
}

/** The tree leaves' order over the members by label; members the tree lacks after them, in data order. */
function treeOrder(model: StaticModel, owner: string, labels: readonly string[]): Int32Array | undefined {
  let newick: string | undefined;
  for (const c of model.spec.children ?? [])
    if (c.kind === 'data')
      for (const n of c.children ?? [])
        if (n.kind === 'tree' && exactAddress(n.params?.of) === owner && typeof n.params?.newick === 'string') {
          newick = n.params.newick as string;
          break;
        }
  if (newick === undefined) return undefined;
  const t = readNewick(newick);
  const byLabel = new Map<string, number>();
  labels.forEach((l, i) => byLabel.has(l) || byLabel.set(l, i));
  const out: number[] = [];
  const seen = new Uint8Array(labels.length);
  for (const leaf of t.leaves) {
    const i = byLabel.get(t.label[leaf]);
    if (i === undefined || seen[i]) continue;
    seen[i] = 1;
    out.push(i);
  }
  for (let i = 0; i < labels.length; i++) if (!seen[i]) out.push(i);
  return Int32Array.from(out);
}

/** A permutation by a field: numbers by value when every value reads as one, text by code points; missing last. */
function fieldOrder(values: readonly unknown[], n: number, descending: boolean): Int32Array {
  const raw: (string | number | undefined)[] = new Array(n);
  let numeric = true;
  for (let i = 0; i < n; i++) {
    const v = values[i];
    raw[i] = v === null || v === undefined || v === '' ? undefined : (v as string | number);
    if (raw[i] !== undefined && !Number.isFinite(Number(raw[i]))) numeric = false;
  }
  if (numeric) for (let i = 0; i < n; i++) if (raw[i] !== undefined) raw[i] = Number(raw[i]);
  const idx = Array.from({ length: n }, (_, i) => i);
  const sign = descending ? -1 : 1;
  idx.sort((a, b) => {
    const x = raw[a];
    const y = raw[b];
    if (x === undefined || y === undefined) return x === y ? a - b : x === undefined ? 1 : -1;
    const c =
      typeof x === 'number' && typeof y === 'number'
        ? x - y
        : String(x) < String(y)
          ? -1
          : String(x) > String(y)
            ? 1
            : 0;
    return sign * c || a - b;
  });
  return Int32Array.from(idx);
}

const Made = new WeakMap<object, Map<string, Members>>();

const keyOf = (o: MemberOrder, n: number) => `${n}|${typeof o === 'string' ? o : `${o.field}|${o.direction}`}`;

/**
 * The members of a categorical axis of `n` members in its order: made once a node, length and order (a reorder is a
 * new order: a new version), kept while the spec's node is.
 */
export function membersOf(model: StaticModel, axis: MemberAxis, n: number): Members {
  const owner = nodeAt(model.spec, axis.owner);
  const holder: object = owner ?? model;
  let byKey = Made.get(holder);
  if (!byKey) Made.set(holder, (byKey = new Map()));
  const key = keyOf(axis.order, n);
  const kept = byKey.get(key);
  if (kept) return kept;
  const labels = labelsOf(owner, axis, n);
  const o = axis.order;
  let order: Int32Array | undefined;
  if (o === 'tree') order = treeOrder(model, axis.owner, labels);
  else if (typeof o === 'object') order = fieldOrder(fieldOf(owner, o.field, labels), n, o.direction === 'descending');
  if (!order) {
    order = new Int32Array(n);
    for (let i = 0; i < n; i++) order[i] = i;
  }
  const made = new Members(order, labels);
  byKey.set(key, made);
  return made;
}
