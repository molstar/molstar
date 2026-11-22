/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { OpenContext, OpenedSource, SourceEntry } from '@molstar/seqx-core/data';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { type LinkSystem, linksBedpe } from './links.js';

// The generated links as a BEDPE file's text, opened by the `bed` runtime as it opens any BEDPE: read whole, indexed by
// both ends, served as the protocol's `pairs` by region (every pair with an end in it) and as bin pairs at its coarser
// levels (the pairs level, made per tile from the held rows). No adapter: what the store reads is the runtime's.

/** A BEDPE source's text, opened by the `bed` runtime. */
export function openBedpeByRegion(
  source: SourceEntry,
  ctx: OpenContext,
  system: LinkSystem,
  text: string = linksBedpe(undefined, undefined, system),
): Promise<OpenedSource> {
  return BedRuntime.open({ ...source, format: 'bedpe', value: text as never }, ctx);
}
