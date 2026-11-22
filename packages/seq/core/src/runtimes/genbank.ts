/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  abortError,
  type Coverage,
  type FormatRuntime,
  type NameHit,
  type OpenContext,
  type OpenedSource,
  type ShapeInfo,
  type SourceEntry,
  type SystemInfo,
  type Tables,
} from '../data/protocol.js';
import { wholeLevel } from '../data/seq-format.js';
import { GenbankCauses } from './genbank-issues.js';
import { parseLocation } from './genbank-location.js';
import { type GenbankRecord, parseGenbank } from './genbank-parse.js';
import { failure, refuse } from './internal/answers.js';
import { lowerBound } from './internal/intervals.js';
import { fieldInfos, type ItemField, type PieceItem, rangesOf, rowsIn, systemItems } from './internal/items.js';
import { decodeText, textBytes } from './internal/text.js';

// GenBank (_spec/README.md "Data", GenBank), the opt-in runtime `@molstar/seqx-core/runtimes/genbank`: a flat file of
// one or more records read whole within `limits` (plain, gzip or BGZF by its magic), each record a sequence named by
// LOCUS's name (its `ACCESSION.VERSION` an alias), `circular` a cyclic system, the molecule its alphabet, ORIGIN's
// letters served as `letters`; its features as `ranges` (the key the `category`, the qualifiers fields, the strand from
// `complement`, `join` and `order` pieces of one item, a join over a circular record's origin one item across it);
// a name index over `/label`, `/gene`, `/locus_tag` and the records' names (EC-228). Positions as the format writes them
// (base 1, ends inclusive), converted at the edge: the source's `base` and `ends` are not read (EC-195).

export { GenbankIssues } from './genbank-issues.js';
export { type Location, MaxDepth, parseLocation } from './genbank-location.js';
export {
  type GenbankFeature,
  type GenbankFile,
  type GenbankLimits,
  type GenbankRecord,
  parseGenbank,
} from './genbank-parse.js';

export const GenbankFormats = ['genbank'] as const;

export interface GenbankOptions {
  /** Records a file (past it the file is refused, `seq.source.too_large`). */
  readonly maxRecords: number;
  /** Qualifiers a feature (past it the feature is left out, said). */
  readonly maxQualifiers: number;
  /** A qualifier's text in characters (past it the feature is left out, said). */
  readonly maxQualifierBytes: number;
  /** Pieces a location (past it the feature is left out, said). */
  readonly maxLocationPieces: number;
}

export const DefaultGenbankOptions: GenbankOptions = {
  maxRecords: 1000,
  maxQualifiers: 256,
  maxQualifierBytes: 64 * 1024,
  maxLocationPieces: 10_000,
};

/** What a host shows of a record besides its sequence: its accession, DEFINITION, organism, molecule. */
export interface GenbankRecordInfo {
  readonly name: string;
  readonly accession?: string;
  readonly definition?: string;
  readonly organism?: string;
  /** LOCUS's molecule type (`DNA`, `mRNA`, `AA`, …). */
  readonly molecule: string;
  readonly length: number;
  readonly circular: boolean;
  readonly features: number;
}

/** An opened GenBank source: its records as a host shows them. */
export interface GenbankSource extends OpenedSource {
  readonly records: readonly GenbankRecordInfo[];
}

/** A feature as an item: its key, its qualifiers, its partial ends. */
export interface GenbankItem extends PieceItem {
  readonly system: string;
  readonly key: string;
  readonly qualifiers: Readonly<Record<string, string>>;
  readonly label: string;
  readonly partialStart: boolean;
  readonly partialEnd: boolean;
}

/** A protein's letters (the language's rule for `protein`). */
const ProteinLetters = /^[A-Za-z*-]*$/;

/** Qualifier names fields at most (the rest said). */
const MaxNames = 256;

/** The qualifiers `findNames` searches, then the records' names. */
const NameQualifiers = ['label', 'gene', 'locus_tag'] as const;
const NameFields = [...NameQualifiers, 'locus', 'accession'];

/** The fields every item has; a qualifier of one of these names is not a field of its own. */
const Fixed = new Set(['label', 'category', 'partial_start', 'partial_end', 'strand']);

const IndexHint = 'A smaller file: split its records';

/** An item's label: `/label`, else `/gene`, else `/locus_tag`, else its key. */
function labelOf(key: string, q: Readonly<Record<string, string>>): string {
  return q.label ?? q.gene ?? q.locus_tag ?? key;
}

/** A record's features as items, each refused one said; ids `<key>@<first>-<last>` (1-based), `#2` on a repeat. */
function itemsOf(r: GenbankRecord, length: number, causes: GenbankCauses, maxPieces: number): GenbankItem[] {
  const out: GenbankItem[] = [];
  const seen = new Map<string, number>();
  const protein = r.alphabet === 'protein';
  for (const f of r.features) {
    const loc = parseLocation(f.location, length, r.circular, maxPieces);
    if (!loc.ok) {
      causes.add('seq.source.genbank.location', `A location not read: ${loc.why}`, `line ${f.line}`);
      continue;
    }
    if (loc.approximate)
      causes.add('seq.source.genbank.approximate', 'A location `a.b` read as its hull', `line ${f.line}`);
    const base = `${f.key}@${loc.start + 1}-${loc.end}`;
    const k = (seen.get(base) ?? 0) + 1;
    seen.set(base, k);
    out.push({
      id: k === 1 ? base : `${base}#${k}`,
      system: r.name,
      start: loc.start,
      end: loc.end,
      strand: protein ? 0 : loc.strand,
      pieces: loc.pieces,
      key: f.key,
      qualifiers: f.qualifiers,
      label: labelOf(f.key, f.qualifiers),
      partialStart: loc.partialStart,
      partialEnd: loc.partialEnd,
    });
  }
  return out;
}

/** The items' fields: the fixed ones, then the qualifiers in order of first use, at most `MaxNames`. */
function itemFields(items: readonly GenbankItem[], causes: GenbankCauses): ItemField<GenbankItem>[] {
  const fixed: ItemField<GenbankItem>[] = [
    { name: 'label', type: 'category', value: (it) => it.label },
    { name: 'category', type: 'category', value: (it) => it.key },
    { name: 'partial_start', type: 'int8', value: (it) => (it.partialStart ? 1 : 0) },
    { name: 'partial_end', type: 'int8', value: (it) => (it.partialEnd ? 1 : 0) },
  ];
  const keys: string[] = [];
  const seen = new Set<string>();
  let past = 0;
  for (const it of items)
    for (const k in it.qualifiers) {
      if (seen.has(k) || Fixed.has(k)) continue;
      seen.add(k);
      if (keys.length < MaxNames) keys.push(k);
      else past++;
    }
  if (past) causes.add('seq.source.genbank.names', `More than ${MaxNames} qualifier names: ${past} not fields`);
  return [
    ...fixed,
    ...keys.map((k) => ({
      name: k,
      type: 'category' as const,
      value: (it: GenbankItem) => it.qualifiers[k] ?? null,
    })),
  ];
}

/** An entry of the name index: a name's item (or a record's whole sequence), the field it was found in. */
interface NameEntry {
  readonly system: string;
  readonly start: number;
  readonly end: number;
  readonly id?: string;
  readonly label: string;
  readonly field: string;
}

/** The name index (EC-228): each lowercased name to its entries, built once at open; exact names first, then prefixes. */
class NameIndex {
  private readonly byName = new Map<string, NameEntry[]>();
  private sorted: string[] | undefined;

  add(name: string | undefined, entry: NameEntry) {
    const key = name?.trim().toLowerCase();
    if (!key) return;
    let list = this.byName.get(key);
    if (!list) this.byName.set(key, (list = []));
    list.push(entry);
  }

  find(text: string, max: number, prefix: boolean): NameHit[] {
    const query = text.trim().toLowerCase();
    const out: NameHit[] = [];
    if (!query || max <= 0) return out;
    const seen = new Set<NameEntry>();
    const add = (entries: readonly NameEntry[]) => {
      for (const e of entries) {
        if (out.length >= max) return;
        if (seen.has(e)) continue;
        seen.add(e);
        out.push({ ...e });
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

async function open(source: SourceEntry, ctx: OpenContext, options: GenbankOptions): Promise<GenbankSource> {
  const at = source.address;
  if (source.format !== 'genbank') refuse(at, 'seq.source.no_data', `Not GenBank: ${source.format}`);
  const causes = new GenbankCauses();
  if (source.base !== undefined || source.ends !== undefined)
    causes.add('seq.source.genbank.transform', 'GenBank fixes base 1 and inclusive ends');

  let bytes: Uint8Array;
  if (typeof source.value === 'string') bytes = new TextEncoder().encode(source.value);
  else if (source.uri !== undefined) {
    try {
      bytes = await ctx.assets.readBytes(source.uri, ctx.signal);
    } catch {
      if (ctx.signal.aborted) throw abortError();
      refuse(at, 'seq.source.unreadable', `The file ${source.uri} could not be read`);
    }
  } else refuse(at, 'seq.source.no_data', 'A GenBank source needs a uri or a text value');
  const text = await textBytes(bytes, ctx.limits.maxBytes, ctx.signal);
  if (!text.ok)
    refuse(at, text.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', text.message, IndexHint);
  if (ctx.signal.aborted) throw abortError();
  const decoded = decodeText(text.bytes);
  if (decoded.replaced) causes.add('seq.source.genbank.encoding', 'Replacement characters in the text');
  const file = parseGenbank(decoded.text, causes, {
    maxRecords: options.maxRecords,
    maxFeatures: ctx.limits.maxRows,
    maxQualifiers: options.maxQualifiers,
    maxQualifierBytes: options.maxQualifierBytes,
    maxLocationPieces: options.maxLocationPieces,
  });
  if (!file.ok)
    refuse(at, file.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed', file.message, IndexHint);

  interface Held {
    readonly record: GenbankRecord;
    readonly length: number;
    readonly items: GenbankItem[];
  }
  const held = new Map<string, Held>();
  const names = new Set<string>();
  for (const r of file.records) {
    if (names.has(r.name)) {
      causes.add('seq.source.genbank.record', `A record named ${r.name} again`, `line ${r.line}`);
      continue;
    }
    if (r.alphabet === 'protein' && r.letters !== undefined && !ProteinLetters.test(r.letters))
      refuse(at, 'seq.source.malformed', `The protein record ${r.name}'s letters are not a protein's`);
    const length = r.letters !== undefined ? r.letters.length : r.length;
    if (r.letters !== undefined && length !== r.length)
      causes.add(
        'seq.source.genbank.length',
        `ORIGIN’s ${length} letters, LOCUS’s length ${r.length}: the letters’`,
        `line ${r.line}`,
      );
    names.add(r.name);
    held.set(r.name, { record: r, length, items: itemsOf(r, length, causes, options.maxLocationPieces) });
  }
  // an accession an alias only where it names no other record
  const aliasOf = (r: GenbankRecord) =>
    r.accession !== undefined && r.accession !== r.name && !names.has(r.accession) ? [r.accession] : undefined;

  const allItems: GenbankItem[] = [];
  for (const h of held.values()) for (const it of h.items) allItems.push(it);
  const fields = itemFields(allItems, causes);
  const systems = new Map([...held].map(([name, h]) => [name, { h, items: systemItems(h.items) }]));

  const index = new NameIndex();
  for (const h of held.values()) {
    const r = h.record;
    const whole = { system: r.name, start: 0, end: h.length, label: r.name };
    index.add(r.name, { ...whole, field: 'locus' });
    if (r.accession !== undefined) index.add(r.accession, { ...whole, field: 'accession' });
    for (const it of h.items)
      for (const q of NameQualifiers)
        index.add(it.qualifiers[q], {
          system: it.system,
          start: it.start,
          end: it.end,
          id: it.id,
          label: it.label,
          field: q,
        });
  }

  let extent = 1;
  const infos: SystemInfo[] = [];
  const shapes: ShapeInfo[] = [];
  const ranges = fieldInfos(fields);
  for (const { h } of systems.values()) {
    const r = h.record;
    const aliases = aliasOf(r);
    extent = Math.max(extent, h.length);
    infos.push({
      name: r.name,
      length: h.length,
      alphabet: r.alphabet,
      ...(r.circular ? { cyclic: true } : {}),
      ...(aliases ? { aliases } : {}),
    });
    if (r.letters !== undefined) shapes.push({ system: r.name, shape: 'letters' });
    shapes.push({ system: r.name, shape: 'ranges', fields: ranges });
  }
  const whole = wholeLevel(extent);
  const coverage: Coverage = { systems: infos, shapes, levels: whole > 1 ? [1, whole] : [1], names: NameFields };
  const records: GenbankRecordInfo[] = [...held.values()].map(({ record: r, length, items }) => ({
    name: r.name,
    ...(r.accession !== undefined ? { accession: r.accession } : {}),
    ...(r.definition !== undefined ? { definition: r.definition } : {}),
    ...(r.organism !== undefined ? { organism: r.organism } : {}),
    molecule: r.molecule,
    length,
    circular: r.circular,
    features: items.length,
  }));
  const maxRows = ctx.limits.maxRows;
  return {
    coverage,
    issues: causes.issues(at),
    records,
    read(request, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      const sys = systems.get(request.system);
      if (!sys) return Promise.resolve(failure('not-served', `No sequence ${request.system}`));
      if (!coverage.levels.includes(request.level))
        return Promise.resolve(failure('not-served', `No level ${request.level}`));
      const tables: { -readonly [K in keyof Tables]: Tables[K] } = {};
      const letters = sys.h.record.letters;
      // letters at level 1 only; a coarse level's tile is its ranges
      if (request.shapes.includes('letters') && letters !== undefined && request.level === 1) {
        const start = Math.max(0, Math.min(request.start, letters.length));
        const end = Math.max(start, Math.min(request.end, letters.length));
        tables.letters = { form: 'letters', start, letters: letters.slice(start, end) };
      }
      if (request.shapes.includes('ranges')) {
        const got = rowsIn(sys.items, request.start, request.end);
        if (got.length / 3 > maxRows)
          return Promise.resolve(failure('too-large', `${got.length / 3} rows in one region`));
        tables.ranges = rangesOf(sys.items, got, fields);
      }
      if (!tables.letters && !tables.ranges)
        return Promise.resolve(failure('not-served', 'The GenBank source serves letters and ranges'));
      return Promise.resolve({ ok: true as const, tables });
    },
    findNames(text, opts, signal) {
      if (signal.aborted) return Promise.reject(abortError());
      return Promise.resolve(index.find(text, Math.max(0, Math.floor(opts.max ?? 20)), opts.prefix ?? false));
    },
    dispose() {},
  };
}

/** A GenBank runtime with its own bounds (`DefaultGenbankOptions` for those not given). */
export function genbankRuntime(over: Partial<GenbankOptions> = {}): FormatRuntime {
  const options = { ...DefaultGenbankOptions, ...over };
  return { name: 'genbank', formats: GenbankFormats, open: (source, ctx) => open(source, ctx, options) };
}

/** The runtime of `genbank` with the default bounds: opt in by listing it. */
export const GenbankRuntime: FormatRuntime = genbankRuntime();
