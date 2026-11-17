/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { SourceIssues } from './issues.js';
import {
  abortError,
  type Coverage,
  type FormatRuntime,
  isAbortError,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type SourceEntry,
} from './protocol.js';
import { SourceOpenError } from './seq-runtime.js';

// FASTA with its index (_spec/README.md "Data", FASTA): a runtime of the protocol. The file is read whole (range reads
// wait for the formats that need them); with an index (`.fai`: name, length, offset, bases a line, bytes a line) the
// records are found by it, after it is checked against the file, so a stale or malformed index is refused, never read
// past. Without one, the index is built from the file. Offsets count bytes: a file that is not ASCII is refused.

/** One line of a `.fai`. */
export interface FaiRecord {
  readonly name: string;
  readonly length: number;
  /** The byte offset of the record's first base. */
  readonly offset: number;
  readonly lineBases: number;
  /** Bytes a line, its line break included. */
  readonly lineWidth: number;
}

const NonAscii = /[^\x00-\x7f]/;

/** A `.fai`'s records, or why it is not one (five fields a line, integers, unique names). */
export function parseFai(text: string): { records?: FaiRecord[]; issue?: string } {
  const records: FaiRecord[] = [];
  const names = new Set<string>();
  const lines = text.split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (!line) continue;
    const f = line.split('\t');
    if (f.length !== 5) return { issue: `line ${i + 1} has ${f.length} fields, not 5` };
    const [name, ...rest] = f;
    const n = rest.map((x) => (/^\d+$/.test(x) ? Number(x) : Number.NaN));
    if (!name || n.some((x) => !Number.isSafeInteger(x)))
      return { issue: `line ${i + 1} is not a name and four integers` };
    if (names.has(name)) return { issue: `the record ${name} is listed twice` };
    names.add(name);
    const [length, offset, lineBases, lineWidth] = n;
    records.push({ name, length, offset, lineBases, lineWidth });
  }
  return { records };
}

/** Where position `p` of a record is in the file. */
function at(r: FaiRecord, p: number): number {
  return r.offset + Math.floor(p / r.lineBases) * r.lineWidth + (p % r.lineBases);
}

/**
 * Why an index does not describe a file, or undefined: each record starts a line after its header, its lines hold
 * `lineBases` bases and a break of `lineWidth − lineBases` bytes, its last base is followed by a line break or the end
 * of the file, no `>` is inside its bytes, and no two records' bytes overlap.
 */
export function checkFai(records: readonly FaiRecord[], fasta: string): string | undefined {
  const spans: [number, number, string][] = [];
  for (const r of records) {
    const why = `the index does not fit the file at ${r.name}`;
    if (r.length > 0 && (r.lineBases < 1 || r.lineWidth <= r.lineBases)) return `${why}: its line sizes`;
    if (r.offset < 1 || r.offset > fasta.length || fasta[r.offset - 1] !== '\n') return `${why}: its offset`;
    const header = fasta.lastIndexOf('>', r.offset - 1);
    const name = header < 0 ? '' : fasta.slice(header + 1, r.offset - 1).split(/\s/)[0];
    if (name !== r.name) return `${why}: its header`;
    if (r.length === 0) continue;
    const last = at(r, r.length - 1);
    if (last >= fasta.length) return `${why}: it ends past the file`;
    const after = fasta.slice(last + 1, last + 3);
    if (after !== '' && after[0] !== '\n' && after !== '\r\n' && !(after === '\r' && last + 2 === fasta.length))
      return `${why}: its last line does not end where it says`;
    if (fasta.slice(r.offset, last + 1).includes('>')) return `${why}: it runs into another record`;
    spans.push([r.offset, last + 1, r.name]);
    // every full line ends where the index says
    for (let p = r.lineBases; p <= r.length; p += r.lineBases) {
      const end = at(r, p - 1) + 1;
      if (p < r.length && !/^\r?\n$/.test(fasta.slice(end, end + r.lineWidth - r.lineBases))) return `${why}: a line`;
    }
  }
  spans.sort((a, b) => a[0] - b[0]);
  for (let i = 1; i < spans.length; i++)
    if (spans[i][0] < spans[i - 1][1]) return `the records ${spans[i - 1][2]} and ${spans[i][2]} overlap`;
  return undefined;
}

/** A record's letters by its index (checked). */
function lettersOf(r: FaiRecord, fasta: string): string {
  let out = '';
  for (let p = 0; p < r.length; p += r.lineBases) {
    const a = at(r, p);
    out += fasta.slice(a, a + Math.min(r.lineBases, r.length - p));
  }
  return out;
}

/**
 * The index of a FASTA file, as `samtools faidx` writes it; throws on lines of unequal length within a record, and on
 * a sequence line after an empty one within a record (an empty line ends its sequence).
 */
export function faidx(fasta: string): string {
  const out: string[] = [];
  const lines = fasta.split('\n');
  let pos = 0;
  let cur:
    | {
        name: string;
        offset: number;
        length: number;
        lineBases: number;
        lineWidth: number;
        short: boolean;
        ended: boolean;
      }
    | undefined;
  const flush = () => {
    if (cur) out.push(`${cur.name}\t${cur.length}\t${cur.offset}\t${cur.lineBases}\t${cur.lineWidth}`);
  };
  for (let i = 0; i < lines.length; i++) {
    const raw = lines[i];
    const width = raw.length + (i < lines.length - 1 ? 1 : 0);
    if (raw.startsWith('>')) {
      flush();
      cur = {
        name: raw.slice(1).split(/\s/)[0],
        offset: pos + width,
        length: 0,
        lineBases: 0,
        lineWidth: 0,
        short: false,
        ended: false,
      };
    } else if (cur && !raw.replace(/\r$/, '').length) {
      if (cur.length > 0) cur.ended = true;
    } else if (cur) {
      if (cur.ended) throw new Error(`A sequence line of ${cur.name} follows an empty line`);
      const bases = raw.replace(/\r$/, '').length;
      if (cur.short || (cur.lineBases && bases > cur.lineBases))
        throw new Error(`The lines of ${cur.name} are of unequal length`);
      if (!cur.lineBases) {
        cur.lineBases = bases;
        cur.lineWidth = width;
      } else if (bases < cur.lineBases) cur.short = true;
      cur.length += bases;
    }
    pos += width;
  }
  flush();
  return out.length ? `${out.join('\n')}\n` : '';
}

function malformed(source: SourceEntry, why: string): SourceOpenError {
  return new SourceOpenError([
    issueOf(SourceIssues, 'seq.source.malformed', `The FASTA source cannot be read: ${why}`, {
      address: source.address,
    }),
  ]);
}

/** A FASTA file and its index (or none: built from the file) as an opened source of letters by record. */
export function openFasta(fasta: string, index?: string, source: SourceEntry = { address: '', format: 'fasta' }) {
  if (NonAscii.test(fasta)) throw malformed(source, 'the file is not ASCII, so byte offsets are not positions');
  let records: readonly FaiRecord[];
  if (index === undefined) {
    try {
      const parsed = parseFai(faidx(fasta));
      if (!parsed.records) throw new Error(parsed.issue);
      records = parsed.records;
    } catch (e) {
      throw malformed(source, String((e as Error).message));
    }
  } else {
    const parsed = parseFai(index);
    if (!parsed.records) throw malformed(source, `the index is malformed: ${parsed.issue}`);
    const why = checkFai(parsed.records, fasta);
    if (why) throw malformed(source, why);
    records = parsed.records;
  }
  const letters = new Map(records.map((r) => [r.name, lettersOf(r, fasta)]));
  const coverage: Coverage = {
    systems: records.map((r) => ({ name: r.name, length: r.length })),
    shapes: records.map((r) => ({ system: r.name, shape: 'letters' as const })),
    levels: [1],
  };
  const opened: OpenedSource = {
    coverage,
    read(request, signal): Promise<RegionAnswer> {
      if (signal.aborted) return Promise.reject(abortError());
      const all = letters.get(request.system);
      if (all === undefined || request.level !== 1 || !request.shapes.includes('letters'))
        return Promise.resolve({
          ok: false,
          kind: 'not-served',
          issue: issueOf(SourceIssues, 'seq.source.not_served', 'The FASTA source serves the letters of its records'),
        });
      const start = Math.min(request.start, all.length);
      const end = Math.min(request.end, all.length);
      return Promise.resolve({
        ok: true,
        tables: { letters: { form: 'letters', start, letters: all.slice(start, end) } },
      });
    },
    dispose() {},
  };
  return opened;
}

async function open(source: SourceEntry, ctx: OpenContext): Promise<OpenedSource> {
  const read = async (uri: string) => {
    try {
      return await ctx.assets.readText(uri, ctx.signal);
    } catch (e) {
      if (isAbortError(e) || ctx.signal.aborted) throw e;
      throw new SourceOpenError([
        issueOf(SourceIssues, 'seq.source.unreadable', `The file ${uri} could not be read`, {
          address: source.address,
        }),
      ]);
    }
  };
  if (source.uri === undefined)
    throw new SourceOpenError([
      issueOf(SourceIssues, 'seq.source.no_data', 'A FASTA source needs a uri', { address: source.address }),
    ]);
  const [fasta, index] = await Promise.all([
    read(source.uri),
    source.index === undefined ? undefined : read(source.index),
  ]);
  return openFasta(fasta, index, source);
}

/** The runtime of FASTA (`format: 'fasta'`, an `index` for its `.fai`). */
export const FastaRuntime: FormatRuntime = { name: 'fasta', formats: ['fasta'], open };
