/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { abortError, type RegionAnswer } from '../data/protocol.js';
import { failure, refuse } from './internal/answers.js';
import { BgzfIssues } from './internal/issues.js';
import { openIndexedLines } from './internal/lines.js';
import { coverageOf, headerSchema, refusedIssues, type VcfOpened, type VcfSource } from './vcf-read.js';
import { readVcfHeader } from './vcf-header.js';
import { VcfCauses } from './vcf-issues.js';
import { parseLine, VcfRowsBuilder, vcfAnswer } from './vcf-parse.js';

// A VCF by region through its tabix index (_spec/README.md "Data", VCF): its header read once at open (io's header
// lines, within `maxHeaderBytes`), each region's lines through the one line parser, at level 1. Loaded when the first
// indexed source opens, so a host reading VCF files whole never loads io's tabix.

export async function openIndexed(o: VcfOpened, uri: string, indexUri: string): Promise<VcfSource> {
  const { source, ctx, options, causes } = o;
  const at = source.address;
  const lines = await openIndexedLines(uri, indexUri, {
    assets: ctx.assets,
    signal: ctx.signal,
    limits: ctx.limits,
    address: at,
    tabix: { maxLineBytes: options.maxLineBytes, maxHeaderBytes: options.maxHeaderBytes, ...options.tabix },
  });
  try {
    const head = await lines.header(ctx.signal);
    // a header past `maxHeaderBytes` is malformed (the format's bound), a range read that failed unreadable
    if (!head.ok)
      refuse(at, head.kind === 'failed' ? 'seq.source.unreadable' : 'seq.source.malformed', `${uri}: ${head.message}`);
    const read = readVcfHeader(head.lines.bytes, options);
    if (!read.ok) refuse(at, 'seq.source.malformed', `${uri}: ${read.message}`);
    const { schema, lengths, refused } = headerSchema(o, read.header);
    const refusedSet = new Set(refused);
    // the index's sequences, and the header's contigs it has no line of (their lengths, no rows)
    const indexed = new Set(lines.index.names);
    const names = [...lengths.keys()];
    for (const name of lines.index.names) if (!lengths.has(name) && !refusedSet.has(name)) names.push(name);
    const issues = [...causes.issues(at), ...refusedIssues(at, refused)];
    if (!lines.eof)
      issues.push(
        issueOf(BgzfIssues, 'seq.source.bgzf.eof', `${uri} has no BGZF EOF block`, {
          address: at,
          hint: BgzfIssues.codes['seq.source.bgzf.eof'].hint,
        }),
      );
    const coverage = coverageOf(names, lengths, schema, [1]);
    const served = new Set(names);
    const maxRows = ctx.limits.maxRows;
    return {
      coverage,
      header: read.header,
      issues,
      async read(request, signal): Promise<RegionAnswer> {
        if (signal.aborted) throw abortError();
        const { system, start, end } = request;
        if (!served.has(system)) return failure('not-served', `No sequence ${system}`);
        if (request.level !== 1) return failure('not-served', 'By region, level 1 only');
        if (!request.shapes.includes('ranges') || !indexed.has(system) || end <= start)
          return vcfAnswer(undefined, schema, request, maxRows);
        const got = await lines.region(system, start, end, signal);
        if (!got.ok) return failure(got.kind, `${uri}, ${system}:${start + 1}-${end}: ${got.message}`);
        const l = got.lines;
        if (l.lineStart.length > maxRows) return failure('too-large', `${l.lineStart.length} lines in one region`);
        // a region's line issues have no place but a failure: the whole read says them
        const quiet = new VcfCauses();
        const rows = new VcfRowsBuilder(schema, lengths.get(system));
        const lineOptions = {
          maxAlleles: options.maxAlleles,
          builder: (s: string) => (s === system ? rows : undefined),
          rowsLeft: () => maxRows - rows.rows,
        };
        const encoding = { replaced: false };
        for (let i = 0; i < l.lineStart.length; i++) {
          const place = `${system}:${l.begin[i] + 1}`;
          if (parseLine(l.bytes, l.lineStart[i], l.lineEnd[i], place, lineOptions, quiet, encoding) < 0)
            return failure('too-large', `More than ${maxRows} alleles in one region`);
        }
        return vcfAnswer(rows.finish(), schema, request, maxRows);
      },
      dispose: () => lines.dispose(),
    };
  } catch (e) {
    lines.dispose();
    throw e;
  }
}
