/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { delayedRuntime, type RegionRequest, SeqRuntime, type StoreClock } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { RegionFile, RegionFormat, RegionSystem } from './region.js';

// The gene's region by region: the region's letters and its transcripts are read from a source that answers region by
// region after a delay; each answer brings the transcripts it overlaps as dynamic nodes, shown a track each by a
// template, and TP53's carry a nested source of their protein's domains, shown on the region through each
// transcript's codon map by a second template (one level down). Nothing about the transcripts is in the spec.

export const title = 'p53’s gene region, by region';

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const source = b.data.source('region', { format: RegionFormat, uri: RegionFile });
  const region = b.data.sequence('tp53-region', { alphabet: 'dna', source, name: RegionSystem });
  b.views.view('main', (v) => {
    const s = v.section('region', { of: region });
    s.track('letters', { label: 'chr17' }).layer('letters', { type: 'letters' });
    s.track('transcripts', {
      label: 'Transcript',
      for_each: { source, node_kind: 'ranges', filter: { map: 'codon' } },
    }).layer('ranges', { type: 'ranges' });
    s.track('domains', {
      label: 'Domains',
      for_each: { source, node_kind: 'ranges', depth: 1 },
    }).layer('ranges', { type: 'ranges' });
  });
  return b.build();
}

export interface RegistryOptions {
  /** The clock the answers wait on: real timers by default. */
  readonly clock?: StoreClock;
  /** Each answer's delay, ms: by default later for regions further along, so they arrive one by one. */
  readonly delayMs?: number | ((request: RegionRequest) => number);
}

/** A region's delay: 300 ms, and 250 more for each 1,024 positions along. */
export const staggered = (r: RegionRequest) => 300 + 250 * Math.floor(r.start / 1024);

/** What the view needs besides the default registry: the delayed `seq` runtime its sources are read with. */
export function registry(options: RegistryOptions = {}): SeqRegistryEntry[] {
  const runtime = delayedRuntime(SeqRuntime, {
    name: 'seq-async',
    formats: [RegionFormat],
    delayMs: options.delayMs ?? staggered,
    ...(options.clock ? { clock: options.clock } : {}),
  });
  return [{ kind: 'format-runtime', runtime }];
}
