/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Column, type ColumnType, categoryColumn } from '@molstar/core/repr/tree';
import type { FieldInfo, RangesField, RangesList } from '../../data/protocol.js';
import { eachMeeting, type IntervalIndex, intervalIndex } from './intervals.js';

// Items of pieces held by sequence, and what a read answers from them (_spec/README.md "Data", Genome runtimes'
// internals): the pieces meeting a region, each with its item's whole extent (`itemStart`, `itemEnd`), and a span row
// (`start = end` at the region's clipped start of the item) for an item meeting the region with no piece in it (a
// region inside an intron), so a tile reads as the whole file clipped. Shared by the GFF and GenBank runtimes, not
// published (`./runtimes/internal/*`).

/** An item of pieces: its identity, extent, strand and pieces (`[start, end, start, end, …]`, sorted and apart). */
export interface PieceItem {
  readonly id: string;
  readonly start: number;
  readonly end: number;
  readonly strand: -1 | 0 | 1;
  readonly pieces: readonly number[];
}

/** A field of the items: its name, type and value of an item (null or NaN missing). */
export interface ItemField<T extends PieceItem = PieceItem> {
  readonly name: string;
  readonly type: ColumnType;
  readonly value: (item: T) => string | number | null;
}

/** One sequence's items: their pieces and the items found by a region (`IntervalIndex` each). */
export interface SystemItems<T extends PieceItem = PieceItem> {
  readonly items: readonly T[];
  readonly pStart: Int32Array;
  readonly pEnd: Int32Array;
  readonly pItem: Int32Array;
  readonly pieces: IntervalIndex;
  readonly iStart: Int32Array;
  readonly iEnd: Int32Array;
  readonly itemsBy: IntervalIndex;
  readonly maxEnd: number;
}

/** A sorted order of `n` values by `key` (stable). */
function orderBy(n: number, key: (i: number) => number): Int32Array {
  const order = Int32Array.from({ length: n }, (_, i) => i);
  order.sort((a, b) => key(a) - key(b) || a - b);
  return order;
}

/** Items of one sequence held for reads. */
export function systemItems<T extends PieceItem>(items: readonly T[]): SystemItems<T> {
  let n = 0;
  for (const it of items) n += it.pieces.length / 2;
  const pStart = new Int32Array(n);
  const pEnd = new Int32Array(n);
  const pItem = new Int32Array(n);
  const iStart = new Int32Array(items.length);
  const iEnd = new Int32Array(items.length);
  let maxEnd = 0;
  for (let i = 0, j = 0; i < items.length; i++) {
    const p = items[i].pieces;
    for (let k = 0; k + 1 < p.length; k += 2, j++) {
      pStart[j] = p[k];
      pEnd[j] = p[k + 1];
      pItem[j] = i;
    }
    iStart[i] = items[i].start;
    iEnd[i] = items[i].end;
    maxEnd = Math.max(maxEnd, items[i].end);
  }
  return {
    items,
    pStart,
    pEnd,
    pItem,
    pieces: intervalIndex(pStart, pEnd),
    iStart,
    iEnd,
    itemsBy: intervalIndex(iStart, iEnd),
    maxEnd,
  };
}

/**
 * The rows of `[start, end)`: each piece meeting it, and a span row for an item meeting it with none, in order of
 * start. `[start, end, item]` triples; empty for an empty region.
 */
export function rowsIn(sys: SystemItems<PieceItem>, start: number, end: number): number[] {
  const out: number[] = [];
  if (end <= start) return out;
  const { pStart, pEnd, pItem, iStart, iEnd } = sys;
  const held = new Set<number>();
  eachMeeting(sys.pieces, pStart, pEnd, start, end, (k) => {
    out.push(pStart[k], pEnd[k], pItem[k]);
    held.add(pItem[k]);
  });
  const spans: number[] = [];
  eachMeeting(sys.itemsBy, iStart, iEnd, start, end, (i) => {
    if (held.has(i)) return;
    const at = Math.max(start, iStart[i]);
    spans.push(at, at, i);
  });
  // in order of start (a long piece or a span row beside the pieces of others)
  if (!spans.length && !sys.pieces.long.length) return out;
  const all = spans.length ? out.concat(spans) : out;
  const order = orderBy(all.length / 3, (r) => all[3 * r]);
  const merged: number[] = [];
  for (const r of order) merged.push(all[3 * r], all[3 * r + 1], all[3 * r + 2]);
  return merged;
}

function column(values: (string | number | null)[], type: ColumnType): Column {
  if (type === 'category') return categoryColumn(values as (string | null)[]);
  const out =
    type === 'int8'
      ? new Int8Array(values.length)
      : type === 'int32'
        ? new Int32Array(values.length)
        : type === 'float32'
          ? new Float32Array(values.length)
          : new Float64Array(values.length);
  for (let r = 0; r < values.length; r++) out[r] = values[r] as number;
  return out as Column;
}

/** The ranges list of rows (`rowsIn`'s triples): identity, strand, fields, every row its item's extent. */
export function rangesOf<T extends PieceItem>(
  sys: SystemItems<T>,
  rows: readonly number[],
  fields: readonly ItemField<T>[],
): RangesList {
  const n = rows.length / 3;
  const start = new Int32Array(n);
  const end = new Int32Array(n);
  const strand = new Int8Array(n);
  const itemStart = new Int32Array(n);
  const itemEnd = new Int32Array(n);
  const ids: string[] = [];
  const values = fields.map(() => [] as (string | number | null)[]);
  for (let r = 0; r < n; r++) {
    const it = sys.items[rows[3 * r + 2]];
    start[r] = rows[3 * r];
    end[r] = rows[3 * r + 1];
    strand[r] = it.strand;
    itemStart[r] = it.start;
    itemEnd[r] = it.end;
    ids.push(it.id);
    for (let f = 0; f < fields.length; f++) values[f].push(fields[f].value(it));
  }
  const out: RangesField[] = fields.map((f, k) => ({ name: f.name, column: column(values[k], f.type) }));
  // the strand a field too, so a layer's data filter and a theme read it (`(item :strand -1)`)
  out.push({ name: 'strand', column: strand });
  return { form: 'list', start, end, strand, id: categoryColumn(ids), fields: out, itemStart, itemEnd };
}

/** The coverage's fields of a list of item fields (and the strand). */
export function fieldInfos(fields: readonly { readonly name: string; readonly type: ColumnType }[]): FieldInfo[] {
  return [...fields.map((f) => ({ name: f.name, type: f.type })), { name: 'strand', type: 'int8' as const }];
}
