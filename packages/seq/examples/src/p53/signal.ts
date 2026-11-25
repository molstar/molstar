/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { BigWigRuntime } from '@molstar/seqx-core/runtimes/bigwig';
import { ConservationFile, ConservationFormat } from './conservation.js';
import { delayedRegistry, type RegistryOptions } from './gene-by-region.js';
import { RegionSystem } from './region.js';

// p53's conservation signal over the gene's region, read from a bigWig: the file's own system and length stand for the
// region (no letters, no other source), and the signal is read by region and zoom level as the view pans and zooms.
// The `bigwig` runtime is opt in, so the view's registry composes it; its answers are delayed as a remote file's are
// (the twin the story and its test share).

export const title = 'p53’s conservation signal from a bigWig';

/** The view's spec. */
export function build(): Node {
  const b = createSeqSpec({ title });
  const source = b.data.source('conservation-bw', {
    format: ConservationFormat,
    uri: ConservationFile,
    attribution: 'A synthetic conservation signal, generated for the examples',
    licence: 'CC0-1.0',
  });
  const region = b.data.sequence('tp53-region', { alphabet: 'dna', source, name: RegionSystem });
  const signal = b.data.values('conservation', { target: region, source });
  b.views.view('main', (v) => {
    const s = v.section('region', { of: region });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('signal', { label: 'Conservation' }).layer('values', { type: 'values', data: signal });
  });
  return b.build();
}

/** What the view needs besides the default registry: the `bigwig` runtime, each answer after `delayMs` (200). */
export function registry(options: RegistryOptions = {}): SeqRegistryEntry[] {
  return delayedRegistry(BigWigRuntime, 'bigwig-async', ConservationFormat, {
    ...options,
    delayMs: options.delayMs ?? 200,
  });
}
