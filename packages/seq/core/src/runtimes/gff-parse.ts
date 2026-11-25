/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GffCauses } from './gff-issues.js';

// GFF3's and GTF's one line parser (_spec/README.md "Data", GFF3 and GTF): nine tab-separated columns into a row, its
// positions converted at the edge (base 1, ends inclusive: fixed by the format, EC-195), its attributes GFF3's
// `key=value;` (percent-decoded, `Parent` and `Alias` lists) or GTF's `key "value";` (unquoted numbers too, a repeated
// key joined by `,`), its identity and parents by the format's own rule (GFF3's `ID` and `Parent`; GTF's `gene_id` and
// `transcript_id`). Directives: `##sequence-region` lengths (EC-188), everything after `##FASTA` ignored but counted.

export type GffFormat = 'gff3' | 'gtf';

/** The largest position a column holds (int32). */
export const MaxPosition = 0x7fffffff;

/** Attributes a row keeps (a file is untrusted input); the rest dropped, said. */
export const MaxAttributes = 256;

/** A line of the file as a row: positions 0-based, end exclusive. */
export interface GffRow {
  readonly seqid: string;
  readonly source: string;
  readonly type: string;
  readonly start: number;
  readonly end: number;
  /** NaN when `.`. */
  readonly score: number;
  readonly strand: -1 | 0 | 1;
  /** 0, 1, 2; −1 when `.`. */
  readonly phase: number;
  /** In file order, values decoded; a GTF key repeated joined by `,`. */
  readonly attributes: Readonly<Record<string, string>>;
  /** GFF3's multi-valued attributes (`MultiValued`) as their values, each decoded after the split. */
  readonly lists?: Readonly<Record<string, readonly string[]>>;
  /** GFF3's `ID`; GTF's `gene_id` on a gene line, `transcript_id` on a transcript line. */
  readonly id?: string;
  /** GFF3's `Parent` list; GTF's `gene_id` of a transcript line, `transcript_id` of the others. */
  readonly parents: readonly string[];
  /** The line's number in the file, 1-based. */
  readonly line: number;
}

/** What a file read whole holds besides its rows. */
export interface GffFile {
  readonly rows: GffRow[];
  /** `##sequence-region` lengths by seqid (EC-188). */
  readonly lengths: Map<string, number>;
  /** Bytes after `##FASTA` (0: none). */
  readonly fastaBytes: number;
  /** Whether a line held U+FFFD. */
  readonly replaced: boolean;
}

/** A GFF3 value: percent-decoded (a malformed escape left as written). */
export function gffDecode(v: string): string {
  if (!v.includes('%')) return v;
  try {
    return decodeURIComponent(v);
  } catch {
    return v;
  }
}

/** A 1-based whole number within int32, else undefined. */
function positionOf(cell: string): number | undefined {
  if (cell.length === 0 || cell.length > 10) return undefined;
  for (let i = 0; i < cell.length; i++) {
    const c = cell.charCodeAt(i);
    if (c < 48 || c > 57) return undefined;
  }
  const v = Number(cell);
  return v >= 1 && v <= MaxPosition ? v : undefined;
}

/** GFF3's attributes of several values (the specification's), split on `,` before each value is decoded. */
export const MultiValued: ReadonlySet<string> = new Set(['Parent', 'Alias', 'Dbxref', 'Ontology_term']);

/** A multi-valued attribute's values: split on `,` (at most `MaxAttributes`), each decoded, empty ones dropped. */
function values(raw: string): string[] {
  const out: string[] = [];
  for (let s = 0; s <= raw.length && out.length < MaxAttributes; ) {
    let e = raw.indexOf(',', s);
    if (e < 0) e = raw.length;
    if (e > s) out.push(gffDecode(raw.slice(s, e)));
    s = e + 1;
  }
  return out;
}

/**
 * GFF3's attributes: `key=value` pairs by `;` (walked, never split whole: a cell of a million `;` allocates no list),
 * decoded; a key written twice keeps the last; the multi-valued ones split on `,` before decoding (an encoded comma,
 * `%2C`, is part of a value). Past `MaxAttributes` keys the walk stops.
 */
function gff3Attributes(
  cell: string,
  keep: (count: number) => boolean,
): { a: Record<string, string>; lists?: Record<string, string[]> } {
  const a: Record<string, string> = Object.create(null);
  let lists: Record<string, string[]> | undefined;
  if (cell === '.' || cell === '') return { a };
  let count = 0;
  for (let s = 0; s < cell.length; ) {
    let e = cell.indexOf(';', s);
    if (e < 0) e = cell.length;
    const eq = cell.indexOf('=', s);
    const at = s;
    s = e + 1;
    if (eq < 0 || eq >= e || eq === at) continue;
    const key = cell.slice(at, eq).trim();
    if (!key) continue;
    if (!(key in a) && !keep(++count)) break;
    const raw = cell.slice(eq + 1, e);
    if (MultiValued.has(key)) {
      const v = values(raw);
      (lists ??= Object.create(null) as Record<string, string[]>)[key] = v;
      a[key] = v.join(',');
    } else a[key] = gffDecode(raw);
  }
  return lists ? { a, lists } : { a };
}

/**
 * GTF's attributes: `key "value";` or `key value;` (a number unquoted), a quoted value read to its closing quote (a `;`
 * inside kept); a key repeated (GENCODE's `tag`) joined by `,`. Undefined when the cell is not attributes.
 */
function gtfAttributes(cell: string, keep: (count: number) => boolean): Record<string, string> | undefined {
  const a: Record<string, string> = Object.create(null);
  if (cell === '.' || cell === '') return a;
  let count = 0;
  const n = cell.length;
  let i = 0;
  const space = (c: number) => c === 32 || c === 9;
  while (i < n) {
    while (i < n && (space(cell.charCodeAt(i)) || cell.charCodeAt(i) === 59)) i++;
    if (i >= n) break;
    const k0 = i;
    while (i < n && !space(cell.charCodeAt(i)) && cell.charCodeAt(i) !== 59) i++;
    const key = cell.slice(k0, i);
    while (i < n && space(cell.charCodeAt(i))) i++;
    let value: string;
    if (cell.charCodeAt(i) === 34) {
      const close = cell.indexOf('"', i + 1);
      if (close < 0) return undefined;
      value = cell.slice(i + 1, close);
      i = close + 1;
    } else {
      const v0 = i;
      while (i < n && cell.charCodeAt(i) !== 59) i++;
      value = cell.slice(v0, i).trim();
      if (!value) return undefined;
    }
    while (i < n && space(cell.charCodeAt(i))) i++;
    if (i < n && cell.charCodeAt(i) !== 59) return undefined;
    if (key in a) a[key] = `${a[key]},${value}`;
    else if (keep(++count)) a[key] = value;
    else break;
  }
  return a;
}

/** Whether a type is a transcript's by GTF's own words (`transcript` and the types `transcriptTypes` lists). */
export type TypeTest = (type: string) => boolean;

/**
 * A file's text read whole into rows. Lines that do not parse are left out (said with their number); a row's attributes
 * past `MaxAttributes` are dropped (said). `maxRows` bounds the data lines: past it, the message why (the caller
 * refuses the source).
 */
export function parseGff(
  text: string,
  format: GffFormat,
  isTranscript: TypeTest,
  causes: GffCauses,
  maxRows: number,
): GffFile | string {
  const rows: GffRow[] = [];
  const lengths = new Map<string, number>();
  const lineCode = format === 'gtf' ? 'seq.source.gtf.line' : 'seq.source.gff.line';
  let fastaBytes = 0;
  let replaced = false;
  const n = text.length;
  let at = 0;
  let line = 0;
  while (at < n) {
    let end = text.indexOf('\n', at);
    if (end < 0) end = n;
    line++;
    let l = text.slice(at, end);
    const next = end + 1;
    if (l.endsWith('\r')) l = l.slice(0, -1);
    if (line === 1 && l.charCodeAt(0) === 0xfeff) l = l.slice(1);
    if (l.charCodeAt(0) === 35) {
      if (l.startsWith('##FASTA')) {
        // its characters, a fair count of its bytes (the bound was checked on the bytes already)
        fastaBytes = n - next;
        break;
      }
      if (l.startsWith('##sequence-region')) {
        const parts = l.slice(17).trim().split(/\s+/);
        const last = parts.length === 3 ? positionOf(parts[2]) : undefined;
        if (last !== undefined) lengths.set(parts[0], last);
      }
      at = next;
      continue;
    }
    at = next;
    if (l.trim() === '') continue;
    if (!replaced && l.includes('�')) {
      replaced = true;
      causes.add('seq.source.gff.encoding', 'Replacement characters (U+FFFD)', `line ${line}`);
    }
    const place = `line ${line}`;
    // the tabs counted before the line is split: a line of a million tabs allocates no list
    let tabs = 0;
    for (let t = l.indexOf('\t'); t >= 0 && tabs < 9; t = l.indexOf('\t', t + 1)) tabs++;
    const cells = tabs === 8 ? l.split('\t') : [];
    if (cells.length !== 9) {
      causes.add(lineCode, `Not nine tab-separated columns`, place);
      continue;
    }
    const s = positionOf(cells[3]);
    const e = positionOf(cells[4]);
    if (s === undefined || e === undefined) {
      causes.add(lineCode, 'A position is not a whole number from 1 to 2^31 − 1', place);
      continue;
    }
    if (e < s) {
      causes.add(lineCode, 'An end before its start', place);
      continue;
    }
    if (rows.length >= maxRows) return `More than ${maxRows} rows`;
    const st = cells[6];
    const strand = st === '+' ? 1 : st === '-' ? -1 : 0;
    const ph = cells[7];
    const phase = ph === '0' ? 0 : ph === '1' ? 1 : ph === '2' ? 2 : -1;
    const score = cells[5] === '.' || cells[5] === '' ? Number.NaN : Number(cells[5]);
    let dropped = false;
    const keep = (count: number) => {
      if (count <= MaxAttributes) return true;
      dropped = true;
      return false;
    };
    const read = format === 'gtf' ? { a: gtfAttributes(cells[8], keep) } : gff3Attributes(cells[8], keep);
    const attributes = read.a;
    const lists = 'lists' in read ? read.lists : undefined;
    if (!attributes) {
      causes.add(lineCode, 'Attributes that are not `key "value";`', place);
      continue;
    }
    if (dropped)
      causes.add('seq.source.gff.attributes', `More than ${MaxAttributes} attributes: the rest dropped`, place);
    const type = format === 'gff3' ? gffDecode(cells[2]) : cells[2];
    let id: string | undefined;
    let parents: string[];
    if (format === 'gff3') {
      id = attributes.ID;
      parents = lists?.Parent ?? [];
    } else if (type === 'gene') {
      id = attributes.gene_id;
      parents = [];
    } else if (isTranscript(type)) {
      id = attributes.transcript_id;
      parents = attributes.gene_id ? [attributes.gene_id] : [];
    } else {
      // a piece's transcript; a line of neither names the gene, so it joins no hierarchy but its transcript's
      parents = attributes.transcript_id ? [attributes.transcript_id] : [];
    }
    rows.push({
      seqid: format === 'gff3' ? gffDecode(cells[0]) : cells[0],
      source: cells[1],
      type,
      start: s - 1,
      end: e,
      score: Number.isNaN(score) ? Number.NaN : score,
      strand,
      phase,
      attributes,
      ...(lists ? { lists } : {}),
      ...(id ? { id } : {}),
      parents,
      line,
    });
  }
  return { rows, lengths, fastaBytes, replaced };
}
