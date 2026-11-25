/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { NameHit } from '../data/protocol.js';
import { type Assembly, type GffItem, mergeIntervals } from './gff-assemble.js';
import type { GffCauses } from './gff-issues.js';
import { type GffRow, MaxAttributes } from './gff-parse.js';
import { lowerBound } from '@molstar/seqx-core/data/intervals';
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

/** The `canonical` field's two values (_spec/README.md "Data", GFF3 and GTF, Tables): what a transcript's tags say. */
export const MANESelect = 'MANE Select';
export const Canonical = 'Canonical';

/**
 * What a transcript's `tag` says of it (GFF3's `tag=a,b`, GTF's repeated `tag "a"` joined by `,`): `MANE Select` for
 * the tag `MANE_Select`, else `Canonical` for `Ensembl_canonical`; null otherwise and for a gene or another feature (a
 * gene item carries its first transcript's attributes: it is no transcript). The viewer never guesses one.
 */
export function canonicalOf(it: GffItem): string | null {
  const tag = it.attributes.tag;
  if (it.role !== 'transcript' || tag === undefined) return null;
  let canonical = false;
  for (const t of tag.split(',')) {
    const name = t.trim();
    if (name === 'MANE_Select') return MANESelect;
    if (name === 'Ensembl_canonical') canonical = true;
  }
  return canonical ? Canonical : null;
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

/** What a source's items are: its transcripts (the default), its genes (transcripts merged) or the file's lines. */
export type GffFeatures = 'transcripts' | 'genes' | 'rows';

/** A whole-sequence feature's types: a feature of one over a whole sequence describes it (Ensembl writes one each). */
export const Describing = new Set(['region', 'chromosome', 'supercontig']);

/**
 * Whether a row describes its sequence (a `region`, `chromosome` or `supercontig` feature from its first position to its
 * length, or to its own end when no length is known, which it then gives): never an item.
 */
export function describes(r: GffRow, lengths: Map<string, number>): boolean {
  if (!Describing.has(r.type) || r.start !== 0) return false;
  const known = lengths.get(r.seqid);
  if (known !== undefined && known !== r.end) return false;
  if (known === undefined) lengths.set(r.seqid, r.end);
  return true;
}

/** Attributes not served as fields: the identity and the hierarchy (`id` and the item's pieces say them). */
const NotFields = new Set(['ID', 'Parent']);

/** The items of the file's lines (`features: rows`): rows sharing an `ID` one item, its pieces those rows, collected
 * once and merged once (linear in the lines). */
export function rowItems(rows: readonly GffRow[]): (GffItem & { readonly row: GffRow })[] {
  const out: (GffItem & { readonly row: GffRow })[] = [];
  /** An `ID`'s item and the pieces its rows gave, merged at the end. */
  const byId = new Map<string, { k: number; start: number; end: number; pieces: number[] }>();
  const seen = new Map<string, number>();
  for (const r of rows) {
    if (r.id !== undefined) {
      const g = byId.get(r.id);
      if (g !== undefined) {
        g.start = Math.min(g.start, r.start);
        g.end = Math.max(g.end, r.end);
        g.pieces.push(r.start, r.end);
        continue;
      }
      byId.set(r.id, { k: out.length, start: r.start, end: r.end, pieces: [r.start, r.end] });
    }
    const base = r.id ?? `${r.type}:${r.start + 1}-${r.end}`;
    const k = (seen.get(base) ?? 0) + 1;
    seen.set(base, k);
    out.push({
      id: k === 1 ? base : `${base}#${k}`,
      seqid: r.seqid,
      type: r.type,
      role: 'feature',
      start: r.start,
      end: r.end,
      strand: r.strand,
      pieces: [r.start, r.end],
      attributes: r.attributes,
      ...(r.lists ? { lists: r.lists } : {}),
      line: r.line,
      row: r,
    });
  }
  for (const g of byId.values())
    if (g.pieces.length > 2) out[g.k] = { ...out[g.k], start: g.start, end: g.end, pieces: mergeIntervals(g.pieces) };
  return out;
}

/** The items a `features` mode serves, from the assembly (or the rows). */
export function itemsOf(features: GffFeatures, a: Assembly, rows: readonly GffRow[]): GffItem[] {
  if (features === 'rows') return rowItems(rows);
  if (features === 'transcripts') return [...a.transcripts, ...a.features];
  // one model a gene: the genes, and the transcripts with none (a gene-level item of their own); no other feature
  const genes = new Set(a.genes.map((g) => g.id));
  return [...a.genes, ...a.transcripts.filter((t) => t.gene === undefined || !genes.has(t.gene))];
}

/** The fields of the items: the fixed ones, then the attributes in order of first use, at most `MaxAttributes`. */
export function itemFields(features: GffFeatures, items: readonly GffItem[], causes: GffCauses): ItemField[] {
  const nan = (v: number | undefined) => (v === undefined ? Number.NaN : v);
  const fixed: ItemField[] = [
    { name: 'label', type: 'category', value: labelOf },
    { name: 'type', type: 'category', value: (it) => it.type },
  ];
  if (features === 'rows') {
    const row = (it: GffItem) => (it as GffItem & { row: GffRow }).row;
    fixed.push(
      { name: 'source', type: 'category', value: (it) => row(it).source },
      { name: 'score', type: 'float32', value: (it) => row(it).score },
      { name: 'phase', type: 'int8', value: (it) => row(it).phase },
    );
  } else
    fixed.push(
      { name: 'gene', type: 'category', value: (it) => it.gene ?? null },
      { name: 'gene_name', type: 'category', value: (it) => it.geneName ?? null },
      { name: 'canonical', type: 'category', value: canonicalOf },
      { name: 'thick_start', type: 'float64', value: (it) => nan(it.thickStart) },
      { name: 'thick_end', type: 'float64', value: (it) => nan(it.thickEnd) },
      { name: 'cds_phase', type: 'int8', value: (it) => it.cdsPhase ?? -1 },
    );
  const taken = new Set([...fixed.map((f) => f.name), 'strand']);
  const keys: string[] = [];
  const seen = new Set<string>();
  let past = 0;
  for (const it of items)
    for (const k in it.attributes) {
      if (seen.has(k) || NotFields.has(k) || taken.has(k)) continue;
      seen.add(k);
      if (keys.length < MaxAttributes) keys.push(k);
      else past++;
    }
  if (past) causes.add('seq.source.gff.attributes', `More than ${MaxAttributes} attribute names: ${past} not fields`);
  return [
    ...fixed,
    ...keys.map((k) => ({ name: k, type: 'category' as const, value: (it: GffItem) => it.attributes[k] ?? null })),
  ];
}
