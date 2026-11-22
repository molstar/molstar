/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import {
  abortError,
  type Coverage,
  type OpenContext,
  type OpenedSource,
  type RegionAnswer,
  type SourceEntry,
} from '../data/protocol.js';
import { assembleTranscripts, type GffItem, PieceTypes, transcriptTest } from './gff-assemble.js';
import { GffCauses } from './gff-issues.js';
import { type GffFormat, type GffRow, parseGff } from './gff-parse.js';
import {
  describes,
  fieldInfos,
  type GffFeatures,
  itemFields,
  itemsOf,
  rangesOf,
  rowsIn,
  systemItems,
} from './gff-read.js';
import type { GffOptions } from './gff.js';
import { failure, refuse } from './internal/answers.js';
import { BgzfIssues } from './internal/issues.js';
import { densityOf } from './internal/levels.js';
import { eachMeeting, intervalIndex } from './internal/intervals.js';
import { indexedLevels, openIndexedLines, type TabixLines } from './internal/lines.js';
import { decodeText } from './internal/text.js';

// GFF3 and GTF by region through a tabix index (_spec/README.md "Data", GFF3 and GTF: by region): a tile's lines
// assembled alone, by type (`assembleTranscripts`), so a tile sees what the whole file does: a transcript's line spans
// its children, so every tile meeting a transcript holds it (and its gene's); a tile inside an intron gives the
// transcript's span row (EC-229). What an item is besides its pieces (its CDS hull, its first CDS's phase, a gene's
// merged transcripts) needs all its children, so a tile is assembled with the lines of the flanks its items reach past
// it (read when the index says they fit with it, and counted in its read cost once a read has found them). Levels from the index's bytes, a density above level 1. Loaded when the first indexed
// source opens, so a host reading GFF files whole never loads io's tabix.

/** An open in progress (`gff.ts`'s). */
export interface GffOpened {
  readonly source: SourceEntry;
  readonly ctx: OpenContext;
  readonly options: GffOptions;
  readonly causes: GffCauses;
  readonly format: GffFormat;
  readonly features: GffFeatures;
}

/** Lines' text, one a line, for the one line parser. */
function textOf(l: TabixLines): string {
  const parts: string[] = [];
  for (let i = 0; i < l.lineStart.length; i++) parts.push(decodeText(l.bytes, l.lineStart[i], l.lineEnd[i]).text);
  return parts.join('\n');
}

/** A flank's lines that do not meet `[start, end)` (those that do were read with the tile), by io's extents. */
function outside(l: TabixLines, start: number, end: number): TabixLines {
  const keep: number[] = [];
  for (let i = 0; i < l.lineStart.length; i++)
    if (!(l.begin[i] < end && Math.max(l.end[i], l.begin[i] + 1) > start)) keep.push(i);
  if (keep.length === l.lineStart.length) return l;
  return {
    bytes: l.bytes,
    lineStart: Uint32Array.from(keep, (i) => l.lineStart[i]),
    lineEnd: Uint32Array.from(keep, (i) => l.lineEnd[i]),
    begin: Int32Array.from(keep, (i) => l.begin[i]),
    end: Int32Array.from(keep, (i) => l.end[i]),
  };
}

const isPiece = new Set(PieceTypes);

/** The piece lines read, found by a region (`intervalIndex`): what an item cut short keeps. */
function pieceLines(rows: readonly GffRow[]) {
  const s: number[] = [];
  const e: number[] = [];
  for (const r of rows)
    if (isPiece.has(r.type)) {
      s.push(r.start);
      e.push(r.end);
    }
  const starts = Int32Array.from(s);
  const ends = Int32Array.from(e);
  const index = intervalIndex(starts, ends);
  return (a: number, b: number) => {
    let met = false;
    eachMeeting(index, starts, ends, a, b, () => {
      met = true;
    });
    return met;
  };
}

/** An item reaching past `[start, end)` with its pieces cut to those meeting a piece line read (none: a span row). */
function piecesRead(it: GffItem, read: (a: number, b: number) => boolean, start: number, end: number): GffItem {
  if (it.role === 'feature' || (it.start >= start && it.end <= end)) return it;
  const kept: number[] = [];
  for (let k = 0; k + 1 < it.pieces.length; k += 2)
    if (read(it.pieces[k], it.pieces[k + 1])) kept.push(it.pieces[k], it.pieces[k + 1]);
  return kept.length === it.pieces.length ? it : { ...it, pieces: kept };
}

export async function openIndexed(o: GffOpened, uri: string, indexUri: string): Promise<OpenedSource> {
  const { source, ctx, options, causes, format, features } = o;
  const at = source.address;
  const lines = await openIndexedLines(uri, indexUri, {
    assets: ctx.assets,
    signal: ctx.signal,
    limits: ctx.limits,
    address: at,
    ...(options.tabix && { tabix: options.tabix }),
  });
  try {
    const isTranscript = transcriptTest(options.transcriptTypes);
    // the header's directives: `##sequence-region` lengths (EC-188)
    const head = await lines.header(ctx.signal);
    if (!head.ok)
      refuse(at, head.kind === 'failed' ? 'seq.source.unreadable' : 'seq.source.malformed', `${uri}: ${head.message}`);
    const parsed = parseGff(textOf(head.lines), format, isTranscript, causes, 0);
    const lengths = typeof parsed === 'string' ? new Map<string, number>() : parsed.lengths;
    const names = [...lengths.keys()];
    for (const name of lines.index.names) if (!lengths.has(name)) names.push(name);
    const indexed = new Set(lines.index.names);
    const issues = causes.issues(at);
    if (!lines.eof)
      issues.push(
        issueOf(BgzfIssues, 'seq.source.bgzf.eof', `${uri} has no BGZF EOF block`, {
          address: at,
          hint: BgzfIssues.codes['seq.source.bgzf.eof'].hint,
        }),
      );
    const levels = indexedLevels(lines, lengths, options);
    // the fixed fields (a tile's attributes are its own, in order of first use there)
    const fixed = fieldInfos(itemFields(features, [], new GffCauses()));
    const coverage: Coverage = {
      systems: names.map((name) => (lengths.has(name) ? { name, length: lengths.get(name)! } : { name })),
      shapes: names.map((system) => ({ system, shape: 'ranges' as const, fields: fixed })),
      levels,
    };
    const served = new Set(names);
    const maxRows = ctx.limits.maxRows;
    const quietRows = (l: readonly TabixLines[], quiet: GffCauses): GffRow[] | string => {
      const file = parseGff(l.map(textOf).join('\n'), format, isTranscript, quiet, maxRows);
      if (typeof file === 'string') return file;
      // a whole-sequence feature describes its sequence, never an item (the whole read's rule, its lengths a copy)
      const known = new Map(lengths);
      return file.rows.filter((r) => !describes(r, known));
    };
    /**
     * The flanks a tile's items reach past it (`[lo, start)`, `[end, hi)`), learned by its last read (a tile's hull is
     * known only once its lines are read): its read cost adds theirs. 512 tiles remembered, the oldest let go.
     */
    const flanks = new Map<string, readonly [number, number]>();
    const remember = (key: string, lo: number, hi: number) => {
      flanks.delete(key);
      flanks.set(key, [lo, hi]);
      if (flanks.size > 512) flanks.delete(flanks.keys().next().value!);
    };
    const flankCost = (system: string, start: number, end: number, lo: number, hi: number) =>
      (lo < start ? lines.cost(system, lo, start) : 0) + (hi > end ? lines.cost(system, end, hi) : 0);
    return {
      coverage,
      issues,
      // the tile's chunks, and those of the flanks its items reach (once a read has found them), from the index
      readCost: (r) => {
        if (!served.has(r.system) || !levels.includes(r.level)) return undefined;
        if (!r.shapes.includes('ranges') || !indexed.has(r.system) || r.end <= r.start) return 0;
        const tile = lines.cost(r.system, r.start, r.end);
        const f = flanks.get(`${r.system}|${r.start}|${r.end}`);
        if (!f) return tile;
        const more = flankCost(r.system, r.start, r.end, f[0], f[1]);
        return tile + more <= lines.maxRegionBytes ? tile + more : tile;
      },
      async read(request, signal): Promise<RegionAnswer> {
        if (signal.aborted) throw abortError();
        const { system, start, end } = request;
        if (!served.has(system)) return failure('not-served', `No sequence ${system}`);
        if (!levels.includes(request.level)) return failure('not-served', `No level ${request.level}`);
        if (!request.shapes.includes('ranges')) return failure('not-served', `The ${format} source serves ranges only`);
        const coarse = request.level > 1 && !request.items;
        if (!indexed.has(system) || end <= start) {
          const empty = rangesOf(systemItems([]), [], itemFields(features, [], new GffCauses()));
          return { ok: true, tables: { ranges: coarse ? densityOf(empty, request) : empty } };
        }
        const place = `${uri}, ${system}:${start + 1}-${end}`;
        const got = await lines.region(system, start, end, signal);
        if (!got.ok) return failure(got.kind, `${place}: ${got.message}`);
        // the hull of the lines meeting the tile (but those describing a sequence): their children are inside it
        let lo = start;
        let hi = end;
        const quiet = new GffCauses();
        let rows = quietRows([got.lines], quiet);
        if (typeof rows === 'string') return failure('too-large', `${place}: ${rows}`);
        for (const r of rows) {
          if (r.seqid !== system) continue;
          if (r.start < lo) lo = r.start;
          if (r.end > hi) hi = r.end;
        }
        // only the flanks are read, when the index says they fit with the tile (`maxRegionBytes`); past it the tile's
        // own lines are assembled (a transcript's line gives its span row, EC-229; its CDS hull, phase and a gene's
        // merged pieces are then the tile's only), so a listed level is never too large
        let partial = false;
        if (lo < start || hi > end) {
          partial = true;
          remember(`${system}|${start}|${end}`, lo, hi);
          if (lines.cost(system, start, end) + flankCost(system, start, end, lo, hi) <= lines.maxRegionBytes) {
            const parts = [got.lines];
            for (const [a, b] of [
              [lo, start],
              [end, hi],
            ] as const) {
              if (a >= b) continue;
              const f = await lines.region(system, a, b, signal);
              if (!f.ok) return failure(f.kind, `${place} (its items to ${a + 1}-${b}): ${f.message}`);
              parts.push(outside(f.lines, start, end));
            }
            rows = quietRows(parts, quiet);
            if (typeof rows === 'string') return failure('too-large', `${place}: ${rows}`);
            partial = false;
          }
        }
        if (rows.length > maxRows) return failure('too-large', `${place}: ${rows.length} lines`);
        const assembly = assembleTranscripts(
          rows.filter((r) => r.seqid === system),
          { transcriptTypes: options.transcriptTypes },
        );
        const all = itemsOf(features, assembly, rows);
        // the flanks not read: an item reaching past the tile keeps only the pieces its read lines give (a transcript of
        // none its span row, never its whole extent as one piece)
        const read = partial ? pieceLines(rows) : undefined;
        const items = read ? all.map((it) => piecesRead(it, read, start, end)) : all;
        const fields = itemFields(features, items, quiet);
        const sys = systemItems(items);
        const triples = rowsIn(sys, start, end);
        if (triples.length / 3 > maxRows) return failure('too-large', `${triples.length / 3} rows in one region`);
        const list = rangesOf(sys, triples, fields);
        return { ok: true, tables: { ranges: coarse ? densityOf(list, request) : list } };
      },
      dispose: () => lines.dispose(),
    };
  } catch (e) {
    lines.dispose();
    throw e;
  }
}
