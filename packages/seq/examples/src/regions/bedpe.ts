/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { issueOf } from '@molstar/core/repr/tree';
import { type OpenContext, type OpenedSource, type SourceEntry, SourceIssues } from '@molstar/seqx-core/data';
import { type BedpeSource, BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { type LinkSystem, linksBedpe } from './links.js';

// The generated links as a BEDPE file's text, opened by the `bed` runtime and served by region as the protocol's `pairs`
// shape (a host adapter's job: the runtime serves a BEDPE's pairs by `pairs(system, start, end)` and lists no shape for
// the store to read): the system's length known, the pairs shape listed at level 1, the store's summaries make the bins.

/** A BEDPE source's text as a source the store reads by region: each region the pairs with an end in it. */
export async function openBedpeByRegion(
  source: SourceEntry,
  ctx: OpenContext,
  system: LinkSystem,
  text: string = linksBedpe(undefined, undefined, system),
): Promise<OpenedSource> {
  const bed = (await BedRuntime.open({ ...source, format: 'bedpe', value: text as never }, ctx)) as BedpeSource;
  return {
    coverage: {
      systems: bed.coverage.systems.map((s) => (s.name === system.name ? { ...s, length: system.length } : s)),
      shapes: bed.coverage.systems.map((s) => ({ system: s.name, shape: 'pairs' as const })),
      levels: [1],
    },
    ...(bed.issues && { issues: bed.issues }),
    async read(request, signal) {
      if (request.level !== 1 || !request.shapes.includes('pairs'))
        return {
          ok: false,
          kind: 'not-served',
          issue: issueOf(SourceIssues, 'seq.source.not_served', 'A BEDPE serves its pairs at level 1'),
        };
      const a = await bed.pairs(request.system, request.start, request.end, signal);
      return a.ok ? { ok: true, tables: { pairs: a.pairs } } : a;
    },
    dispose() {
      bed.dispose();
    },
  };
}
