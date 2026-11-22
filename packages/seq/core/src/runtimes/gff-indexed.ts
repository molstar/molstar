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
import { assembleTranscripts, transcriptTest } from './gff-assemble.js';
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
import { indexedLevels, openIndexedLines, type TabixLines } from './internal/lines.js';
import { decodeText } from './internal/text.js';

// GFF3 and GTF by region through a tabix index (_spec/README.md "Data", GFF3 and GTF: by region): a tile's lines
// assembled alone, by type (`assembleTranscripts`), so a tile sees what the whole file does: a transcript's line spans
// its children, so every tile meeting a transcript holds it (and its gene's); a tile inside an intron gives the
// transcript's span row (EC-229). What an item is besides its pieces (its CDS hull, its first CDS's phase, a gene's
// merged transcripts) needs all its children, so a tile is assembled from the hull of the lines meeting it (one more
// read when they reach past it). Levels from the index's bytes, a density above level 1. Loaded when the first indexed
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
    const quietRows = (l: TabixLines, quiet: GffCauses): GffRow[] | string => {
      const file = parseGff(textOf(l), format, isTranscript, quiet, maxRows);
      if (typeof file === 'string') return file;
      // a whole-sequence feature describes its sequence, never an item (the whole read's rule, its lengths a copy)
      const known = new Map(lengths);
      return file.rows.filter((r) => !describes(r, known));
    };
    return {
      coverage,
      issues,
      readCost: (r) =>
        served.has(r.system) && levels.includes(r.level)
          ? r.shapes.includes('ranges') && indexed.has(r.system) && r.end > r.start
            ? lines.cost(r.system, r.start, r.end)
            : 0
          : undefined,
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
        let got = await lines.region(system, start, end, signal);
        if (!got.ok) return failure(got.kind, `${place}: ${got.message}`);
        // the hull of the lines meeting the tile (but those describing a sequence): their children are inside it
        let lo = start;
        let hi = end;
        const quiet = new GffCauses();
        let rows = quietRows(got.lines, quiet);
        if (typeof rows === 'string') return failure('too-large', `${place}: ${rows}`);
        for (const r of rows) {
          if (r.seqid !== system) continue;
          if (r.start < lo) lo = r.start;
          if (r.end > hi) hi = r.end;
        }
        if (lo < start || hi > end) {
          got = await lines.region(system, lo, hi, signal);
          if (!got.ok) return failure(got.kind, `${place} (its items to ${lo + 1}-${hi}): ${got.message}`);
          rows = quietRows(got.lines, quiet);
          if (typeof rows === 'string') return failure('too-large', `${place}: ${rows}`);
        }
        if (rows.length > maxRows) return failure('too-large', `${place}: ${rows.length} lines`);
        const assembly = assembleTranscripts(
          rows.filter((r) => r.seqid === system),
          { transcriptTypes: options.transcriptTypes },
        );
        const items = itemsOf(features, assembly, rows);
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
