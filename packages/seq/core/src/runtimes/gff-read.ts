/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Column, type ColumnType, categoryColumn } from '@molstar/core/repr/tree';
import type { FieldInfo, NameHit, RangesField, RangesList } from '../data/protocol.js';
import type { GffItem } from './gff-assemble.js';

// A GFF3's or GTF's items held by sequence, and what a read answers from them (_spec/README.md "Data", GFF3 and GTF):
// the pieces meeting a region, each with its item's whole extent (`itemStart`, `itemEnd`), and a span row (`start =
// end` at the region's clipped start of the item) for an item meeting the region with no piece in it (a region inside
// an intron), so a tile reads as the whole file clipped; and the name index (EC-228), built once at open.

/** A field of the items: its name, type and value of an item (null or NaN missing). */
export interface ItemField {
  readonly name: string;
  readonly type: ColumnType;
  readonly value: (item: GffItem) => string | number | null;
}

/** One sequence's items: pieces sorted by start, each with its item; items sorted by start. */
export interface SystemItems {
  readonly items: readonly GffItem[];
  readonly pStart: Int32Array;
  readonly pEnd: Int32Array;
  readonly pItem: Int32Array;
  readonly maxPiece: number;
  /** Items by start, and the longest item (an item meeting `start` starts at `start − maxItem` or later). */
  readonly iOrder: Int32Array;
  readonly maxItem: number;
  readonly maxEnd: number;
}

/** A sorted order of `n` values by `key` (stable). */
function orderBy(n: number, key: (i: number) => number): Int32Array {
  const order = Int32Array.from({ length: n }, (_, i) => i);
  order.sort((a, b) => key(a) - key(b) || a - b);
  return order;
}

/** Items of one sequence held for reads. */
export function systemItems(items: readonly GffItem[]): SystemItems {
  let n = 0;
  for (const it of items) n += it.pieces.length / 2;
  const s: number[] = [];
  const e: number[] = [];
  const owner: number[] = [];
  let maxPiece = 1;
  let maxItem = 1;
  let maxEnd = 0;
  for (let i = 0; i < items.length; i++) {
    const p = items[i].pieces;
    for (let k = 0; k + 1 < p.length; k += 2) {
      s.push(p[k]);
      e.push(p[k + 1]);
      owner.push(i);
      maxPiece = Math.max(maxPiece, p[k + 1] - p[k]);
    }
    maxItem = Math.max(maxItem, items[i].end - items[i].start);
    maxEnd = Math.max(maxEnd, items[i].end);
  }
  const order = orderBy(n, (i) => s[i]);
  const pStart = new Int32Array(n);
  const pEnd = new Int32Array(n);
  const pItem = new Int32Array(n);
  for (let k = 0; k < n; k++) {
    pStart[k] = s[order[k]];
    pEnd[k] = e[order[k]];
    pItem[k] = owner[order[k]];
  }
  const iOrder = orderBy(items.length, (i) => items[i].start);
  return { items, pStart, pEnd, pItem, maxPiece, iOrder, maxItem, maxEnd };
}

/** The first index of `sorted` (by `at`) whose value is at least `v`. */
function lowerBound(count: number, at: (k: number) => number, v: number): number {
  let lo = 0;
  let hi = count;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (at(mid) < v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

const meets = (s: number, e: number, start: number, end: number) => s < end && Math.max(e, s + 1) > start;

/**
 * The rows of `[start, end)`: each piece meeting it, and a span row for an item meeting it with none, in order of
 * start. `[start, end, item]` triples; empty for an empty region.
 */
export function rowsIn(sys: SystemItems, start: number, end: number): number[] {
  const out: number[] = [];
  if (end <= start) return out;
  const { pStart, pEnd, pItem, items, iOrder } = sys;
  const held = new Set<number>();
  for (let k = lowerBound(pStart.length, (k) => pStart[k], start - sys.maxPiece); k < pStart.length; k++) {
    if (pStart[k] >= end) break;
    if (!meets(pStart[k], pEnd[k], start, end)) continue;
    out.push(pStart[k], pEnd[k], pItem[k]);
    held.add(pItem[k]);
  }
  const spans: number[] = [];
  for (let k = lowerBound(iOrder.length, (k) => items[iOrder[k]].start, start - sys.maxItem); k < iOrder.length; k++) {
    const i = iOrder[k];
    const it = items[i];
    if (it.start >= end) break;
    if (held.has(i) || !meets(it.start, it.end, start, end)) continue;
    const at = Math.max(start, it.start);
    spans.push(at, at, i);
  }
  if (!spans.length) return out;
  // merged by start (a span row beside the pieces of others)
  const all = [...out, ...spans];
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
export function rangesOf(sys: SystemItems, rows: readonly number[], fields: readonly ItemField[]): RangesList {
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
export function fieldInfos(fields: readonly ItemField[]): FieldInfo[] {
  return [...fields.map((f) => ({ name: f.name, type: f.type })), { name: 'strand', type: 'int8' as const }];
}

/** An entry of the name index: a name's item, the field it was found in. */
interface NameEntry {
  readonly item: GffItem;
  readonly field: string;
  readonly value: string;
}

/**
 * The name index (EC-228): every value of `fields` (a list value split at `,`) of every item, lowercased, to its
 * items; built once. `find` gives the exact names first, then those starting with the text when asked, at most `max`,
 * each item once a field; a hit's `id` the item's when `served(item)` (an item of the source's `features`), else none
 * (its range addresses it).
 */
export class NameIndex {
  private readonly byName = new Map<string, NameEntry[]>();
  private sorted: string[] | undefined;

  constructor(
    items: Iterable<GffItem>,
    readonly fields: readonly string[],
    private readonly served: (item: GffItem) => boolean,
  ) {
    for (const item of items)
      for (const field of fields) {
        const v = item.attributes[field] ?? (field === 'ID' && item.role !== 'feature' ? item.id : undefined);
        if (v === undefined) continue;
        for (const name of field === 'Alias' ? v.split(',') : [v]) {
          const key = name.trim().toLowerCase();
          if (!key) continue;
          let list = this.byName.get(key);
          if (!list) this.byName.set(key, (list = []));
          if (!list.some((e) => e.item === item && e.field === field)) list.push({ item, field, value: name.trim() });
        }
      }
  }

  get size() {
    return this.byName.size;
  }

  find(text: string, max: number, prefix: boolean): NameHit[] {
    const query = text.trim().toLowerCase();
    const out: NameHit[] = [];
    if (!query || max <= 0) return out;
    const seen = new Set<string>();
    const add = (entries: readonly NameEntry[]) => {
      for (const e of entries) {
        if (out.length >= max) return;
        const it = e.item;
        const key = `${it.id}|${it.seqid}|${it.start}|${it.end}|${e.field}`;
        if (seen.has(key)) continue;
        seen.add(key);
        const label = it.attributes.Name ?? it.attributes.transcript_name ?? it.attributes.gene_name ?? e.value;
        out.push({
          system: it.seqid,
          start: it.start,
          end: it.end,
          ...(this.served(it) ? { id: it.id } : {}),
          label,
          field: e.field,
        });
      }
    };
    const exact = this.byName.get(query);
    if (exact) add(exact);
    if (!prefix || out.length >= max) return out;
    const keys = (this.sorted ??= [...this.byName.keys()].sort());
    for (let k = lowerBound(keys.length, (k) => (keys[k] < query ? -1 : 1), 0); k < keys.length; k++) {
      const key = keys[k];
      if (!key.startsWith(query)) break;
      if (key === query) continue;
      add(this.byName.get(key)!);
      if (out.length >= max) break;
    }
    return out;
  }
}
