/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue, issueOf } from '@molstar/core/repr/tree';
import { abortError, type Coverage, type OpenedSource, type RegionAnswer } from '../data/protocol.js';
import { BedCauses } from './bed-issues.js';
import { cellsOf, fieldInfos } from './bed-parse.js';
import { answerOf, bedgraphIntervals, bedItems, dataCells, type Opened, schemaOf, spanLines } from './bed-read.js';
import { failure, refuse } from './internal/answers.js';
import { BgzfIssues } from './internal/issues.js';
import { indexedLevels, openIndexedLines } from './internal/lines.js';
import { decodeText } from './internal/text.js';

// A BED or bedGraph by region through its tabix index (_spec/README.md "Data", BED family): each region's lines through
// the one line parser, at the levels the index's bytes allow (a BED's density, a bedGraph's bins above level 1, items
// when asked), each tile's read cost the index's. Loaded when the first indexed source opens, so a host reading BED files whole never
// loads io's BGZF reader and tabix.

export async function openIndexed(o: Opened, uri: string, indexUri: string): Promise<OpenedSource> {
  const { source, ctx, causes } = o;
  const format = o.format as 'bed' | 'bedgraph';
  const lines = await openIndexedLines(uri, indexUri, {
    assets: ctx.assets,
    signal: ctx.signal,
    limits: ctx.limits,
    address: source.address,
    tabix: o.options.tabix,
  });
  try {
    const issues: Issue[] = causes.issues(source.address);
    if (!lines.eof)
      issues.push(
        issueOf(BgzfIssues, 'seq.source.bgzf.eof', `${uri} has no BGZF EOF block`, {
          address: source.address,
          hint: BgzfIssues.codes['seq.source.bgzf.eof'].hint,
        }),
      );
    const names = lines.index.names;
    // a BED's columns are its first line's (BED `N`): read at open, so the coverage lists its fields
    let fields: Coverage['shapes'][number]['fields'];
    let schema: ReturnType<typeof schemaOf> | undefined;
    if (format === 'bed') {
      const head = await lines.head(ctx.signal);
      if (!head.ok)
        refuse(
          source.address,
          head.kind === 'too-large' ? 'seq.source.too_large' : 'seq.source.malformed',
          `${uri}: ${head.message}`,
        );
      schema = schemaOf(source, head.line ? cellsOf(decodeText(head.line).text).cells.length : 3);
      fields = fieldInfos(schema);
    }
    const shape = format === 'bed' ? ('ranges' as const) : ('values' as const);
    const levels = indexedLevels(lines, new Map(), o.options);
    const coverage: Coverage = {
      systems: names.map((name) => ({ name })),
      shapes: names.map((system) => (fields ? { system, shape, fields } : { system, shape })),
      levels,
    };
    const maxRows = ctx.limits.maxRows;
    return {
      coverage,
      issues,
      // a tile's read cost: the chunks io reads for it, from the position before (below), nothing read
      readCost: (request) =>
        names.includes(request.system) && levels.includes(request.level)
          ? request.shapes.includes(shape) && request.end > request.start
            ? lines.cost(request.system, Math.max(0, request.start - 1), request.end)
            : 0
          : undefined,
      async read(request, signal): Promise<RegionAnswer> {
        if (signal.aborted) throw abortError();
        const { system, start, end } = request;
        if (!names.includes(system)) return failure('not-served', `No sequence ${system}`);
        if (!levels.includes(request.level)) return failure('not-served', `No level ${request.level}`);
        if (!request.shapes.includes(shape)) return answerOf(format, undefined, request, maxRows, true);
        // from a position before: io meets a line where it overlaps (an item of no length at `start` would be missed), the
        // parser's rule then keeps exactly what the whole read does (a zero-length item at its start)
        const got = await lines.region(system, Math.max(0, start - 1), end, signal);
        if (!got.ok) return failure(got.kind, `${uri}, ${system}:${start}-${end}: ${got.message}`);
        const l = got.lines;
        if (l.lineStart.length > maxRows) return failure('too-large', `${l.lineStart.length} lines in one region`);
        // a region's line issues have no place but a failure: the whole read says them
        const quiet = new BedCauses();
        const rows = dataCells(
          spanLines(l.bytes, l.lineStart, l.lineEnd, (i) => `${system}:${l.begin[i] + 1}`, quiet),
          quiet,
        );
        if (schema) {
          const items = bedItems(rows, schema, quiet, maxRows);
          if (typeof items === 'string') return failure('too-large', items);
          return answerOf(format, items.get(system), request, maxRows, true);
        }
        return answerOf(format, bedgraphIntervals(rows, quiet).get(system), request, maxRows, true);
      },
      dispose: () => lines.dispose(),
    };
  } catch (e) {
    lines.dispose();
    throw e;
  }
}
