/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { NameHit } from '../data/protocol.js';
import type { GffItem } from './gff-assemble.js';
import { lowerBound } from './internal/intervals.js';
import type { ItemField as GenericField, SystemItems as GenericItems } from './internal/items.js';

export { fieldInfos, rangesOf, rowsIn, systemItems } from './internal/items.js';

// A GFF3's or GTF's items held by sequence, and what a read answers from them (_spec/README.md "Data", GFF3 and GTF):
// the items and their answers are the shared ones (`internal/items.ts`: pieces with their item's extent, span rows), so
// a tile reads as the whole file clipped; and the name index (EC-228), built once at open.

/** A field of GFF items. */
export type ItemField = GenericField<GffItem>;
/** One sequence's GFF items. */
export type SystemItems = GenericItems<GffItem>;

/** An item's label: GFF3's `Name`, a gene's `gene_name` or a transcript's `transcript_name` (GTF), else its id. */
export function labelOf(it: GffItem): string {
  const a = it.attributes;
  return a.Name ?? (it.role === 'gene' ? a.gene_name : (a.transcript_name ?? a.gene_name)) ?? it.id;
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
        // GTF writes its gene's names on every line: a transcript's `gene_name` and `gene_id` name its gene, not it
        if ((field === 'gene_name' || field === 'gene_id') && item.role === 'transcript' && item.gene !== undefined)
          continue;
        const v = item.attributes[field] ?? (field === 'ID' && item.role !== 'feature' ? item.id : undefined);
        if (v === undefined) continue;
        // a multi-valued attribute by its values (split before decoding: an encoded comma is part of a name)
        for (const name of item.lists?.[field] ?? [v]) {
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
        const label = labelOf(it);
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
