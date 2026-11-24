/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { GenbankCauses } from './genbank-issues.js';

// A GenBank flat file's records (_spec/README.md "Data", GenBank): `LOCUS` … `//`, each its name, length, molecule and
// topology, DEFINITION, ACCESSION and VERSION, SOURCE's ORGANISM, its feature table (a key, a location's text, the
// qualifiers) and ORIGIN's letters (digits and spaces dropped, case kept). The other header fields skipped, said once.
// One pass over the lines; every scan bounded to its line, every text it holds to its bound before it is held; never a
// throw for a file's content.

/** The bounds of a read (a file is untrusted input). */
export interface GenbankLimits {
  /** Records a file. */
  readonly maxRecords: number;
  /** Features a file (the store's `limits.maxRows`). */
  readonly maxFeatures: number;
  /** Qualifiers a feature. */
  readonly maxQualifiers: number;
  /** A qualifier's text (a repeated one's together), and a DEFINITION's, in characters. */
  readonly maxQualifierBytes: number;
  /** Pieces a location; its text at most 64 characters a piece. */
  readonly maxLocationPieces: number;
}

/** A sequence's letters, at most (int32 positions). */
export const MaxLetters = 0x7fffffff;

/** A feature as written: its key, its location's text (whitespace dropped), its qualifiers (a repeat joined by `; `). */
export interface GenbankFeature {
  readonly key: string;
  readonly location: string;
  readonly qualifiers: Readonly<Record<string, string>>;
  /** A repeated qualifier's values each (its `qualifiers` value them joined). */
  readonly lists?: Readonly<Record<string, readonly string[]>>;
  /** Its first line, 1-based. */
  readonly line: number;
}

/** A record as written. */
export interface GenbankRecord {
  /** LOCUS's name. */
  readonly name: string;
  /** LOCUS's length. */
  readonly length: number;
  /** LOCUS's molecule type (`DNA`, `RNA`, `mRNA`, `AA` for a protein's `aa`). */
  readonly molecule: string;
  readonly alphabet: 'dna' | 'rna' | 'protein';
  readonly circular: boolean;
  /** VERSION's `ACCESSION.VERSION`, else ACCESSION's first. */
  readonly accession?: string;
  readonly definition?: string;
  readonly organism?: string;
  /** ORIGIN's letters as read (absent: no ORIGIN). */
  readonly letters?: string;
  readonly features: readonly GenbankFeature[];
  /** Its LOCUS line, 1-based. */
  readonly line: number;
}

export type GenbankFile =
  | { readonly ok: true; readonly records: readonly GenbankRecord[] }
  | { readonly ok: false; readonly kind: 'too-large' | 'malformed'; readonly message: string };

/** Header keywords read; the others skipped (said once a keyword). */
const Read = new Set(['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'SOURCE', 'FEATURES', 'ORIGIN', 'BASE', '//']);
/** Keywords of an assembly of other records: the record not read. */
const Assembly = new Set(['CONTIG', 'WGS', 'WGS_SCAFLD', 'TSA', 'TLS']);
/** A feature key (INSDC's: letters, digits, `_`, `-`, `'`, `*`). */
const FeatureKey = /^[A-Za-z0-9_'*-]+$/;
/** A location's first characters: a position, a partial end, an operator. */
const LocationStart = /^(?:[<>]?\d|complement\(|join\(|order\()/;

/** A feature being read. */
interface OpenFeature {
  key: string;
  location: string[];
  locationBytes: number;
  /** Each name's values (a repeat a value more) and their characters together, joined once at its end. */
  qualifiers: Map<string, { values: string[]; bytes: number }>;
  count: number;
  line: number;
  /** Why it is left out, once known. */
  bad?: { code: 'seq.source.genbank.feature' | 'seq.source.genbank.location'; why: string };
}

/** A qualifier being read. */
interface OpenQualifier {
  name: string;
  parts: string[];
  bytes: number;
  /** Inside its quotes, past this line. */
  open: boolean;
  /** Not held (`/translation`): only its quotes followed. */
  skip: boolean;
  /** Its line, 1-based. */
  line: number;
}

/** A record being read. */
interface OpenRecord {
  name: string;
  length: number;
  molecule: string;
  alphabet: 'dna' | 'rna' | 'protein';
  circular: boolean;
  accession?: string;
  version?: string;
  definition?: string;
  organism?: string;
  letters?: string[];
  letterCount: number;
  features: GenbankFeature[];
  line: number;
  /** Why it is not read, once known. */
  refused?: string;
}

/** Where a line is read: before a record, its header, its feature table, its letters. */
type State = 'outside' | 'header' | 'features' | 'origin';

const isSpace = (c: number) => c === 32 || c === 9;

/**
 * Where a quoted value's closing quote is in `s` from `from` (`""` an escaped quote), else −1: a scan of the text it
 * is given, nothing past it.
 */
function closingQuote(s: string, from: number): number {
  for (let k = from; k < s.length; k++) {
    if (s.charCodeAt(k) !== 34) continue;
    if (s.charCodeAt(k + 1) === 34) k++;
    else return k;
  }
  return -1;
}

/** The records of a GenBank flat file, or why the file is refused (a bound passed: `too-large`). */
export function parseGenbank(text: string, causes: GenbankCauses, limits: GenbankLimits): GenbankFile {
  const records: GenbankRecord[] = [];
  let featureCount = 0;
  let sawLocus = false;
  let rec: OpenRecord | undefined;
  let state: State = 'outside';
  /** The header keyword a continuation line belongs to. */
  let keyword = '';
  let feature: OpenFeature | undefined;
  let qualifier: OpenQualifier | undefined;
  const maxLocation = limits.maxLocationPieces * 64;
  let tooLarge: string | undefined;

  const closeQualifier = () => {
    const q = qualifier;
    qualifier = undefined;
    if (!q || !feature || feature.bad) return;
    // a quote left open: its text to here is the value (said), the feature kept
    if (q.open)
      causes.add('seq.source.genbank.feature', 'A quote left open: closed where its value ends', `line ${q.line}`);
    if (q.skip) return;
    let v = q.parts.join(' ');
    if (v.startsWith('"')) {
      const end = closingQuote(v, 1);
      v = (end < 0 ? v.slice(1) : v.slice(1, end)).replaceAll('""', '"');
    }
    let at = feature.qualifiers.get(q.name);
    if (!at) feature.qualifiers.set(q.name, (at = { values: [], bytes: 0 }));
    at.bytes += v.length + (at.values.length ? 2 : 0);
    if (at.bytes > limits.maxQualifierBytes) {
      feature.bad = {
        code: 'seq.source.genbank.feature',
        why: `A qualifier past ${limits.maxQualifierBytes} characters`,
      };
      return;
    }
    at.values.push(v);
  };
  const closeFeature = () => {
    closeQualifier();
    const f = feature;
    feature = undefined;
    if (!f || !rec) return;
    if (f.bad) {
      causes.add(f.bad.code, f.bad.why, `line ${f.line}`);
      return;
    }
    if (++featureCount > limits.maxFeatures) {
      tooLarge = `More than ${limits.maxFeatures} features`;
      return;
    }
    // a record without a prototype: a qualifier named `__proto__` is a name like any other
    const qualifiers: Record<string, string> = Object.create(null);
    let lists: Record<string, string[]> | undefined;
    for (const [name, q] of f.qualifiers) {
      qualifiers[name] = q.values.join('; ');
      if (q.values.length > 1) (lists ??= Object.create(null) as Record<string, string[]>)[name] = q.values;
    }
    rec.features.push({
      key: f.key,
      location: f.location.join(''),
      qualifiers,
      ...(lists ? { lists } : {}),
      line: f.line,
    });
  };
  const closeRecord = () => {
    closeFeature();
    const r = rec;
    rec = undefined;
    state = 'outside';
    if (!r) return;
    if (r.refused !== undefined) {
      causes.add('seq.source.genbank.record', r.refused, `line ${r.line}`);
      return;
    }
    if (records.length >= limits.maxRecords) {
      tooLarge = `More than ${limits.maxRecords} records`;
      return;
    }
    const accession = r.version ?? r.accession;
    records.push({
      name: r.name,
      length: r.length,
      molecule: r.molecule,
      alphabet: r.alphabet,
      circular: r.circular,
      ...(accession !== undefined ? { accession } : {}),
      ...(r.definition !== undefined ? { definition: r.definition } : {}),
      ...(r.organism !== undefined ? { organism: r.organism } : {}),
      ...(r.letters ? { letters: r.letters.join('') } : {}),
      features: r.features,
      line: r.line,
    });
  };

  /** The first whitespace-separated token of `[s, e)`, and where it ends. */
  const token = (s: number, e: number): [string, number] => {
    while (s < e && isSpace(text.charCodeAt(s))) s++;
    let k = s;
    while (k < e && !isSpace(text.charCodeAt(k))) k++;
    return [text.slice(s, k), k];
  };

  const locus = (s: number, e: number, line: number): OpenRecord => {
    const words = text.slice(s, e).split(/[ \t]+/);
    const r: OpenRecord = {
      name: words[1] ?? '',
      length: 0,
      molecule: 'DNA',
      alphabet: 'dna',
      circular: false,
      letterCount: 0,
      features: [],
      line,
    };
    const unit = words.findIndex((w, k) => k > 1 && (w === 'bp' || w === 'aa'));
    const length = unit > 1 && /^\d{1,10}$/.test(words[unit - 1]) ? Number(words[unit - 1]) : Number.NaN;
    if (!r.name || unit === 2 || !(length >= 0 && length <= MaxLetters)) {
      r.refused = 'A LOCUS line without a name or a length';
      return r;
    }
    r.length = length;
    if (words[unit] === 'aa') {
      r.molecule = 'AA';
      r.alphabet = 'protein';
    } else {
      const m = words[unit + 1];
      if (m !== undefined && m !== 'linear' && m !== 'circular') r.molecule = m;
      r.alphabet = /RNA/i.test(r.molecule) ? 'rna' : 'dna';
    }
    for (let k = unit + 1; k < words.length; k++) if (words[k].toLowerCase() === 'circular') r.circular = true;
    return r;
  };

  /** A header line (no indent: a keyword; an indent: a sub-keyword or a continuation). */
  const header = (s: number, e: number, line: number, indent: number) => {
    if (!rec) return;
    if (indent === 0) {
      const [word, after] = token(s, e);
      keyword = word;
      if (word === '//') return closeRecord();
      if (word === 'FEATURES') {
        state = 'features';
        return;
      }
      if (word === 'ORIGIN') {
        state = 'origin';
        rec.letters = [];
        return;
      }
      if (Assembly.has(word)) {
        rec.refused ??= `A ${word} assembly of other records`;
        return;
      }
      if (!Read.has(word)) {
        causes.add('seq.source.genbank.header', 'Header fields not read', word.slice(0, 32));
        return;
      }
      const rest = text.slice(after, e).trim();
      if (word === 'DEFINITION') rec.definition = rest.slice(0, limits.maxQualifierBytes);
      else if (word === 'ACCESSION') rec.accession = token(after, e)[0] || undefined;
      else if (word === 'VERSION') rec.version = token(after, e)[0] || undefined;
      return;
    }
    const [word, after] = token(s, e);
    if (indent < 12 && word === 'ORGANISM' && keyword === 'SOURCE') {
      rec.organism = text.slice(after, e).trim().slice(0, limits.maxQualifierBytes);
      return;
    }
    // a DEFINITION at its bound takes no more lines (each would copy it again)
    if (
      indent >= 12 &&
      keyword === 'DEFINITION' &&
      rec.definition !== undefined &&
      rec.definition.length < limits.maxQualifierBytes
    ) {
      const more = `${rec.definition} ${text.slice(s, e).trim()}`;
      rec.definition = more.slice(0, limits.maxQualifierBytes);
    }
  };

  /** A line of the feature table (indented). */
  const featureLine = (s: number, e: number, line: number, indent: number) => {
    const c = text.charCodeAt(s + indent);
    // inside a quoted value every line is its text, a `/` at its start too, but a line at the key's indent with a key
    // and a location: the quote was left open (closed there, said), the line a feature of its own
    const keyLine = qualifier?.open === true && indent < 21 && c !== 47 && isFeatureStart(s, e);
    if ((qualifier?.open && !keyLine) || (indent >= 21 && c !== 47) || (c === 47 && feature)) {
      if (!feature) {
        causes.add('seq.source.genbank.feature', 'A feature table line before any feature', `line ${line}`);
        return;
      }
      const content = text.slice(s + indent, e).trimEnd();
      if (qualifier?.open || (qualifier && c !== 47)) return qualifierText(content);
      if (c === 47) return startQualifier(content, line);
      // the location continued: whitespace dropped
      feature.locationBytes += content.length;
      if (feature.locationBytes > maxLocation)
        feature.bad ??= {
          code: 'seq.source.genbank.location',
          why: `A location of more than ${maxLocation} characters`,
        };
      if (!feature.bad) feature.location.push(content.replace(/\s+/g, ''));
      return;
    }
    if (c === 47) {
      causes.add('seq.source.genbank.feature', 'A qualifier before any feature', `line ${line}`);
      return;
    }
    closeFeature();
    const [key, after] = token(s, e);
    const loc = text.slice(after, e).trim();
    feature = { key: key.slice(0, 64), location: [], locationBytes: loc.length, qualifiers: new Map(), count: 0, line };
    if (loc.length > maxLocation)
      feature.bad = { code: 'seq.source.genbank.location', why: `A location of more than ${maxLocation} characters` };
    else feature.location.push(loc.replace(/\s+/g, ''));
  };

  /** Whether `[s, e)` is a feature's first line: a key, then a location's start (a position, `<`, `>`, an operator). */
  const isFeatureStart = (s: number, e: number): boolean => {
    const [key, after] = token(s, e);
    if (!key || key.length > 64 || !FeatureKey.test(key)) return false;
    const [loc] = token(after, Math.min(e, after + 64));
    return LocationStart.test(loc);
  };

  const startQualifier = (content: string, line: number) => {
    closeQualifier();
    const f = feature!;
    // `/name=value` or `/name`: the `=` looked for in this line only
    let eq = -1;
    for (let k = 1; k < content.length; k++)
      if (content.charCodeAt(k) === 61) {
        eq = k;
        break;
      }
    const name = (eq < 0 ? content.slice(1) : content.slice(1, eq)).trim().slice(0, 64);
    const value = eq < 0 ? 'true' : content.slice(eq + 1);
    const skip = name === 'translation';
    if (skip) causes.add('seq.source.genbank.translation', '`/translation` not read', `line ${line}`);
    else if (++f.count > limits.maxQualifiers)
      f.bad ??= { code: 'seq.source.genbank.feature', why: `More than ${limits.maxQualifiers} qualifiers` };
    const open = value.startsWith('"') && closingQuote(value, 1) < 0;
    qualifier = { name, parts: [], bytes: 0, open, skip, line };
    hold(value);
  };

  const qualifierText = (content: string) => {
    const q = qualifier!;
    if (q.open) q.open = closingQuote(content, 0) < 0;
    hold(content.trim());
  };

  const hold = (part: string) => {
    const q = qualifier!;
    if (q.skip || feature!.bad) return;
    q.bytes += part.length + 1;
    if (q.bytes > limits.maxQualifierBytes + 2) {
      feature!.bad = {
        code: 'seq.source.genbank.feature',
        why: `A qualifier past ${limits.maxQualifierBytes} characters`,
      };
      return;
    }
    q.parts.push(part);
  };

  const originLine = (s: number, e: number) => {
    const r = rec!;
    const part = text.slice(s, e).replace(/[\s\d]+/g, '');
    if (!part) return;
    r.letterCount += part.length;
    if (r.letterCount > MaxLetters) {
      tooLarge = `A record of more than ${MaxLetters} letters`;
      return;
    }
    r.letters!.push(part);
  };

  let lineNo = 0;
  for (let pos = 0; pos < text.length && tooLarge === undefined; ) {
    lineNo++;
    const nl = text.indexOf('\n', pos);
    const next = nl < 0 ? text.length : nl + 1;
    let e = nl < 0 ? text.length : nl;
    if (e > pos && text.charCodeAt(e - 1) === 13) e--;
    const s = pos;
    pos = next;
    let indent = 0;
    while (s + indent < e && text.charCodeAt(s + indent) === 32) indent++;
    if (s + indent === e) continue;
    if (indent === 0 && text.startsWith('LOCUS', s) && isSpace(text.charCodeAt(s + 5))) {
      sawLocus = true;
      if (rec) {
        rec.refused ??= 'A record cut short (a LOCUS before its `//`)';
        closeRecord();
      }
      rec = locus(s, e, lineNo);
      state = 'header';
      keyword = 'LOCUS';
      continue;
    }
    // the state the closures above set (a cast: the compiler sees no assignment through them)
    const now = state as State;
    if (now === 'outside') continue;
    if (now === 'origin') {
      if (indent === 0 && text.startsWith('//', s)) closeRecord();
      else if (indent === 0) {
        // a keyword after the letters (CONTIG after ORIGIN): a header line again
        state = 'header';
        header(s, e, lineNo, 0);
      } else originLine(s, e);
      continue;
    }
    // a line from the first column ends the table, a quote left open too (its feature then left out)
    if (now === 'features' && indent > 0) {
      featureLine(s, e, lineNo, indent);
      continue;
    }
    if (now === 'features') {
      closeFeature();
      state = 'header';
    }
    header(s, e, lineNo, indent);
  }
  if (tooLarge !== undefined) return { ok: false, kind: 'too-large', message: tooLarge };
  if (!sawLocus) return { ok: false, kind: 'malformed', message: 'No LOCUS line: not a GenBank file' };
  if (rec) {
    // the file ends inside a record: the records before it kept
    const cut = rec;
    rec = undefined;
    causes.add(
      'seq.source.genbank.record',
      'The last record cut short (no `//` at the end of the file)',
      `line ${cut.line}`,
    );
  }
  if (!records.length) return { ok: false, kind: 'malformed', message: 'No complete record (`LOCUS` … `//`)' };
  return { ok: true, records };
}
