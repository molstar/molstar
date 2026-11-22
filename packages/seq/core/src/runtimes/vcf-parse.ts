/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { type Column, type ColumnType, categoryColumn } from '@molstar/core/repr/tree';
import type { FieldInfo, RangesField, RangesList, RegionAnswer, RegionRequest } from '../data/protocol.js';
import { failure } from './internal/answers.js';
import { eachMeeting, type IntervalIndex, intervalIndex } from './internal/intervals.js';
import { decodeText } from './internal/text.js';
import type { VcfHeader } from './vcf-header.js';
import type { VcfCauses } from './vcf-issues.js';

// The VCF's one line parser (_spec/README.md "Data", VCF): a data line's eight fixed columns (the samples' never
// decoded) into a row an ALT allele, each with the line's extent by the index's rule (`POS − 1` to `INFO/END` when INFO
// has `END=` past `POS − 1`, else to `POS − 1 + len(REF)`; htslib's `tbx_parse1`, io's `readRegion`), the anchor base
// kept; read whole or from a tabix region alike, so a region read equals the whole read clipped.

/** The largest position (int32). */
export const MaxPosition = 0x7fffffff;

/** An `Integer` INFO value that is missing: BCF's missing int32 (`0x80000000`). */
export const MissingInt = -0x80000000;

/** A field's declared type (the source's `types`); `string` is held as a category. */
type DeclaredType = 'int32' | 'float32' | 'float64' | 'category' | 'string';

/** The fixed fields every row has, in order: the label and `name` the VCF `ID`, then the alleles and the record's. */
export const FixedFields: readonly { readonly name: string; readonly type: ColumnType }[] = [
  { name: 'label', type: 'category' },
  { name: 'name', type: 'category' },
  { name: 'ref', type: 'category' },
  { name: 'alt', type: 'category' },
  { name: 'qual', type: 'float32' },
  { name: 'filter', type: 'category' },
  { name: 'category', type: 'category' },
];

/** An INFO key read as a field: under `name` (the key, or a `fields` mapping's), by its `Number` and type. */
export interface InfoField {
  readonly key: string;
  readonly name: string;
  /** `A` (an allele's value), `R` (the reference's and the allele's), `0` (a flag), `1`, or another (its text). */
  readonly number: string;
  readonly type: ColumnType;
}

export interface VcfSchema {
  readonly info: readonly InfoField[];
  /** An INFO key's field. */
  readonly byKey: ReadonlyMap<string, number>;
  /** The keys the header declares or the source names, fields or not (a key outside them is said once). */
  readonly known: ReadonlySet<string>;
}

/**
 * The fields of a header's INFO keys: those the source's `fields` (a field's name to its key) or `types` (a key's type)
 * name first, then the header's order, at most `maxInfoFields`, the rest said once with their count; a key named as a
 * fixed field said once and not read. `Integer` int32 and `Float` float32 when one value an allele (`Number=1` or `A`),
 * `Flag` and every other a category (`Number=R` the reference's and the allele's value, other numbers their text);
 * `types` overrides a field's type. A key `fields` or `types` names that the header does not declare is read as text.
 */
export function vcfSchema(
  header: VcfHeader,
  fields: Readonly<Record<string, string>> | undefined,
  types: Readonly<Record<string, DeclaredType>> | undefined,
  maxInfoFields: number,
  causes: VcfCauses,
): VcfSchema {
  const declared = new Map(header.info.map((d) => [d.id, d]));
  const nameOf = new Map<string, string>();
  for (const [name, key] of Object.entries(fields ?? {})) if (!nameOf.has(key)) nameOf.set(key, name);
  const keys: string[] = [];
  const seen = new Set<string>();
  const want = (key: string) => {
    if (!seen.has(key)) {
      seen.add(key);
      keys.push(key);
    }
  };
  for (const key of nameOf.keys()) want(key);
  for (const t of Object.keys(types ?? {})) if (declared.has(t) || !FixedFields.some((f) => f.name === t)) want(t);
  for (const d of header.info) want(d.id);
  const fixed = new Set(FixedFields.map((f) => f.name));
  const info: InfoField[] = [];
  let past = 0;
  for (const key of keys) {
    const name = nameOf.get(key) ?? key;
    if (fixed.has(name)) {
      causes.add('seq.source.vcf.cell', `An INFO key named as a fixed field, not read`, key);
      continue;
    }
    if (info.length >= maxInfoFields) {
      past++;
      continue;
    }
    const d = declared.get(key);
    const number = d ? (d.type === 'Flag' ? '0' : d.number) : '.';
    const single = number === '1' || number === 'A';
    const own: ColumnType =
      d && single && d.type === 'Integer' ? 'int32' : d && single && d.type === 'Float' ? 'float32' : 'category';
    const t = types?.[name] ?? types?.[key];
    const type: ColumnType = t === undefined ? own : t === 'string' ? 'category' : t;
    info.push({ key, name, number, type });
  }
  if (past) causes.add('seq.source.vcf.info', `INFO keys past the ${maxInfoFields} read: ${past} not fields`);
  return { info, byKey: new Map(info.map((f, k) => [f.key, k])), known: seen };
}

/** The coverage's fields of a schema's tables. */
export function fieldInfos(schema: VcfSchema): FieldInfo[] {
  return [...FixedFields, ...schema.info.map((f) => ({ name: f.name, type: f.type }))];
}

/** The values of the `category` field, in the order a key lists them (a categorical palette's `values`). */
export const VcfCategories = ['SNV', 'MNV', 'insertion', 'deletion', 'indel', 'structural', 'breakend'] as const;

/** A variant allele's category by its text and the reference's (EC-196). */
export function alleleCategory(ref: string, alt: string): string {
  if (alt.startsWith('<')) return 'structural';
  if (alt.includes('[') || alt.includes(']') || alt.startsWith('.') || alt.endsWith('.')) return 'breakend';
  if (alt.length === ref.length) return alt.length === 1 ? 'SNV' : 'MNV';
  if (alt.length > ref.length && alt.startsWith(ref)) return 'insertion';
  if (alt.length < ref.length && ref.startsWith(alt)) return 'deletion';
  return 'indel';
}

/** An allele as an id's part: itself, or `<length>:<FNV-1a hex>` past 32 characters (a sequence-resolved SV's). */
export function alleleKey(allele: string): string {
  if (allele.length <= 32) return allele;
  const codes = new Uint8Array(allele.length);
  for (let i = 0; i < allele.length; i++) codes[i] = allele.charCodeAt(i) & 0xff;
  return `${allele.length}:${hashFnv32a(codes).toString(16).padStart(8, '0')}`;
}

/** An integer at the start of `text` from `at`, as `strtoll` reads it (a sign, digits): NaN without digits. */
function leadingInteger(text: string, at: number): number {
  let i = at;
  let sign = 1;
  if (text[i] === '+' || text[i] === '-') {
    if (text[i] === '-') sign = -1;
    i++;
  }
  const d0 = i;
  let v = 0;
  while (i < text.length && text.charCodeAt(i) >= 48 && text.charCodeAt(i) <= 57 && v <= Number.MAX_SAFE_INTEGER / 10)
    v = v * 10 + text.charCodeAt(i++) - 48;
  return i === d0 ? Number.NaN : sign * v;
}

/**
 * A line's extent by the index's rule (io's, htslib's `tbx_parse1`): `POS − 1` to `INFO/END` when INFO starts with `END=`
 * or holds `;END=` with a number past `POS − 1` (whatever the alleles), else to `POS − 1 + len(REF)`.
 */
export function vcfExtent(pos: number, ref: string, info: string): [number, number] {
  const begin = Math.max(0, pos - 1);
  let end = begin + ref.length;
  const at = info.startsWith('END=') ? 4 : ((i) => (i < 0 ? -1 : i + 5))(info.indexOf(';END='));
  if (at >= 0 && at < info.length && info[at] !== '.') {
    const v = leadingInteger(info, at);
    if (v > begin) end = v;
  }
  return [begin, end];
}

/** Whether an ALT allele is no variant: missing, an upstream deletion, gVCF's reference block. */
function noVariant(alt: string): string | undefined {
  if (alt === '.' || alt === '') return 'An ALT `.`';
  if (alt === '*') return 'An ALT `*` (an upstream deletion)';
  if (alt === '<NON_REF>' || alt === '<*>') return `An ALT ${alt} (a reference block)`;
  return undefined;
}

/** A cell as its field's type: a number (NaN missing, or `MissingInt`), else its text (null missing). */
function cellValue(text: string | undefined, type: ColumnType, place: string, causes: VcfCauses) {
  const missing = text === undefined || text === '' || text === '.';
  if (type === 'category') return missing ? null : text;
  if (
    type === 'int32' ||
    type === 'int8' ||
    type === 'int16' ||
    type === 'uint8' ||
    type === 'uint16' ||
    type === 'uint32'
  ) {
    if (missing) return MissingInt;
    const v = /^[+-]?\d+$/.test(text) ? Number(text) : Number.NaN;
    if (Number.isInteger(v) && v > MissingInt && v <= MaxPosition) return v;
    causes.add('seq.source.vcf.cell', 'An integer field’s cell is no int32: read as missing', place);
    return MissingInt;
  }
  if (missing) return Number.NaN;
  const v = Number(text);
  if (Number.isNaN(v) && !/^nan$/i.test(text))
    causes.add('seq.source.vcf.cell', 'A number field’s cell is no number: read as missing', place);
  return v;
}

/** A system's rows as a read builds them: a row an allele, in file order. */
export class VcfRowsBuilder {
  readonly start: number[] = [];
  readonly end: number[] = [];
  readonly key: string[] = [];
  /** The fixed fields' values, then the INFO fields'. */
  readonly values: (string | number | null)[][];
  private readonly seen = new Map<string, number>();

  constructor(
    readonly schema: VcfSchema,
    /** The system's length when known: rows past it left out or clipped (EC-51). */
    readonly length: number | undefined,
  ) {
    this.values = Array.from({ length: FixedFields.length + schema.info.length }, () => []);
  }

  get rows() {
    return this.start.length;
  }

  finish(): VcfRows {
    const start = Int32Array.from(this.start);
    const end = Int32Array.from(this.end);
    let maxEnd = 0;
    for (let i = 0; i < end.length; i++) if (end[i] > maxEnd) maxEnd = end[i];
    return { start, end, key: this.key, values: this.values, index: intervalIndex(start, end), maxEnd };
  }

  /** A row of an allele of a parsed line. */
  push(begin: number, end: number, key: string, values: readonly (string | number | null)[]) {
    const n = (this.seen.get(key) ?? 0) + 1;
    this.seen.set(key, n);
    this.start.push(begin);
    this.end.push(end);
    this.key.push(n === 1 ? key : `${key}#${n}`);
    for (let f = 0; f < values.length; f++) this.values[f].push(values[f]);
  }
}

/** A system's rows held for reads: by start (`index`), each with its id and its fields' values. */
export interface VcfRows {
  readonly start: Int32Array;
  readonly end: Int32Array;
  readonly key: readonly string[];
  readonly values: readonly (readonly (string | number | null)[])[];
  readonly index: IntervalIndex;
  readonly maxEnd: number;
}

export interface LineOptions {
  readonly maxAlleles: number;
  /** The builder of a line's system, or undefined when its system is not read (another sequence of a region). */
  readonly builder: (system: string) => VcfRowsBuilder | undefined;
  /** Rows the read may still add (`limits.maxRows` less those made): past it the line is not read, false returned. */
  readonly rowsLeft: () => number;
}

const tab = 9;

/** Where a line's eighth column ends (its INFO's end), or −1 when it has fewer than eight columns. */
function fixedEnd(bytes: Uint8Array, s: number, e: number): number {
  let at = s;
  for (let k = 0; k < 7; k++) {
    const t = bytes.indexOf(tab, at);
    if (t < 0 || t >= e) return -1;
    at = t + 1;
  }
  const t = bytes.indexOf(tab, at);
  return t < 0 || t > e ? e : t;
}

/**
 * Parses a data line `bytes[s, e)` (its eight fixed columns decoded, the samples never): its alleles' rows into its
 * system's builder. A line that does not parse is left out (`.line`), one with more ALT alleles than `maxAlleles` too
 * (`.alleles`, counted from its commas before they are split); an allele that is no variant makes no row (`.allele`).
 * The rows made, or −1 when the rows left would be passed (the read is `too-large`).
 */
export function parseLine(
  bytes: Uint8Array,
  s: number,
  e: number,
  place: string,
  options: LineOptions,
  causes: VcfCauses,
  encoding: { replaced: boolean },
): number {
  const stop = fixedEnd(bytes, s, e);
  if (stop < 0) {
    causes.add('seq.source.vcf.line', 'Fewer than eight columns', place);
    return 0;
  }
  const d = decodeText(bytes, s, stop);
  if (d.replaced && !encoding.replaced) {
    encoding.replaced = true;
    causes.add('seq.source.vcf.encoding', 'Replacement characters (U+FFFD)', place);
  }
  const [chrom, posText, id, ref, altText, qualText, filterText, infoText] = d.text.split('\t');
  const builder = options.builder(chrom);
  if (!builder) return 0;
  if (!/^\d+$/.test(posText) || Number(posText) > MaxPosition) {
    causes.add('seq.source.vcf.line', 'A position that is not a whole number to 2^31 − 1', place);
    return 0;
  }
  if (!ref || ref === '.') {
    causes.add('seq.source.vcf.line', 'No REF', place);
    return 0;
  }
  const pos = Number(posText);
  let [begin, end] = vcfExtent(pos, ref, infoText);
  if (end > MaxPosition) {
    causes.add('seq.source.vcf.line', 'An end past 2^31 − 1', place);
    return 0;
  }
  const schema = builder.schema;
  let commas = 0;
  for (let i = altText.indexOf(','); i >= 0; i = altText.indexOf(',', i + 1)) commas++;
  if (commas + 1 > options.maxAlleles) {
    causes.add('seq.source.vcf.alleles', `More than ${options.maxAlleles} ALT alleles`, place);
    return 0;
  }
  const length = builder.length;
  if (length !== undefined && (begin >= length || end > length)) {
    causes.add('seq.source.vcf.range', 'Past the contig’s length: left out, or clipped', place);
    if (begin >= length) return 0;
    end = length;
  }
  const alts = altText.split(',');
  let variants = 0;
  for (const alt of alts) if (!noVariant(alt)) variants++;
  if (variants > options.rowsLeft()) return -1;

  // INFO's values of the schema's keys, as text (a flag's: '')
  const info: (string | undefined)[] = new Array(schema.info.length);
  if (infoText && infoText !== '.') {
    for (let i = 0; i < infoText.length; ) {
      let j = infoText.indexOf(';', i);
      if (j < 0) j = infoText.length;
      const eq = infoText.indexOf('=', i);
      const key = eq >= 0 && eq < j ? infoText.slice(i, eq) : infoText.slice(i, j);
      const f = schema.byKey.get(key);
      if (f !== undefined) info[f] ??= eq >= 0 && eq < j ? infoText.slice(eq + 1, j) : '';
      else if (key && key !== 'END' && !schema.known.has(key))
        causes.add('seq.source.vcf.info', 'An INFO key the header does not declare', key);
      i = j + 1;
    }
  }
  const name = id === '.' || id === '' ? null : id;
  const qual = qualText === '.' || qualText === undefined ? Number.NaN : Number(qualText);
  if (Number.isNaN(qual) && qualText !== '.' && qualText !== undefined)
    causes.add('seq.source.vcf.cell', 'A QUAL that is no number: read as missing', place);
  const filter = filterText === '.' || filterText === '' ? null : filterText;
  const values: (string | number | null)[] = new Array(FixedFields.length + schema.info.length);
  const refKey = alleleKey(ref);
  let made = 0;
  for (let k = 0; k < alts.length; k++) {
    const alt = alts[k];
    const why = noVariant(alt);
    if (why) {
      causes.add('seq.source.vcf.allele', why, place);
      continue;
    }
    values[0] = name;
    values[1] = name;
    values[2] = ref;
    values[3] = alt;
    values[4] = qual;
    values[5] = filter;
    values[6] = alleleCategory(ref, alt);
    for (let f = 0; f < schema.info.length; f++) {
      const field = schema.info[f];
      const text = info[f];
      let cell: string | undefined;
      if (field.number === '0') cell = text === undefined ? undefined : 'true';
      else if (text === undefined) cell = undefined;
      else if (field.number === 'A') cell = text.split(',')[k];
      else if (field.number === 'R') {
        const v = text.split(',');
        cell = v[0] === undefined || v[k + 1] === undefined ? undefined : `${v[0]},${v[k + 1]}`;
      } else cell = text;
      values[FixedFields.length + f] = cellValue(cell, field.type, place, causes);
    }
    builder.push(begin, end, `${pos}:${refKey}>${alleleKey(alt)}`, values);
    made++;
  }
  return made;
}

/** The rows of `[start, end)` (each meeting it), in order of start, then of the file. */
export function rowsIn(rows: VcfRows, start: number, end: number): number[] {
  const out: number[] = [];
  if (end <= start) return out;
  eachMeeting(rows.index, rows.start, rows.end, start, end, (i) => out.push(i));
  return out.sort((a, b) => rows.start[a] - rows.start[b] || a - b);
}

function column(values: readonly (string | number | null)[], type: ColumnType, rows: readonly number[]): Column {
  if (type === 'category') return categoryColumn(rows.map((r) => values[r] as string | null));
  const out =
    type === 'int32'
      ? new Int32Array(rows.length)
      : type === 'float64'
        ? new Float64Array(rows.length)
        : new Float32Array(rows.length);
  for (let i = 0; i < rows.length; i++) out[i] = values[rows[i]] as number;
  return out as Column;
}

/** The answer of a system's rows over a request (whole and region reads alike): `ranges` lists at any level read. */
export function vcfAnswer(
  rows: VcfRows | undefined,
  schema: VcfSchema,
  request: RegionRequest,
  maxRows: number,
): RegionAnswer {
  if (!request.shapes.includes('ranges')) return failure('not-served', 'A VCF source serves ranges only');
  const got = rows ? rowsIn(rows, request.start, request.end) : [];
  if (got.length > maxRows) return failure('too-large', `${got.length} rows in one region`);
  const names = [...FixedFields, ...schema.info];
  const fields: RangesField[] = names.map((f, k) => ({
    name: f.name,
    column: rows ? column(rows.values[k], f.type, got) : column([], f.type, []),
  }));
  const list: RangesList = {
    form: 'list',
    start: rows ? Int32Array.from(got, (r) => rows.start[r]) : new Int32Array(0),
    end: rows ? Int32Array.from(got, (r) => rows.end[r]) : new Int32Array(0),
    id: categoryColumn(rows ? got.map((r) => rows.key[r]) : []),
    fields,
  };
  return { ok: true, tables: { ranges: list } };
}
