/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { delayedRuntime, type StoreClock } from '@molstar/seqx-core/data';
import { createSeqSpec, ideogramNavigator } from '@molstar/seqx-core/lang';
import { BedFormats, BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { withRecordedPages } from '../recorded-pages.js';
import {
  Chr17Length,
  CytobandsFile,
  GeneratedBedFile,
  GeneratedLength,
  GeneratedSystem,
  MisreadBedFile,
} from './files.js';

// The BED family's views (Features/Genome formats/BED): each a file read by the opt-in `bed` runtime, its system's
// length from an inline chrom.sizes (the same runtime). hg38's cytobands (the recorded real file) as chr17's bands and
// its ideogram; generated BED12 genes, a block an exon; a misread file and its issues.

/** The sections' addresses (what the tests read). */
export const CytobandsSection = '/views/main/chr17';
export const GenesSection = '/views/main/genes';
export const MisreadSection = '/views/main/misread';

/** chr17's bands from UCSC's cytobands: `stain` and `label` by their column numbers (`bed 4+1`). */
export function cytobands(): Node {
  const b = createSeqSpec({ title: 'hg38 cytobands from BED' });
  const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: `chr17\t${Chr17Length}\n` as never });
  const bed = b.data.source('cytobands', {
    format: 'bed',
    uri: CytobandsFile,
    fields: { stain: '4', label: '3' },
    attribution: 'UCSC Genome Browser, hg38 cytoBandIdeo',
    licence: 'LicenseRef-UCSC-freely-usable',
  });
  const chr = b.data.sequence('chr17', { alphabet: 'dna', source: sizes, name: 'chr17' });
  const bands = b.data.ranges('bands', { target: chr, source: bed });
  b.views.view('main', (v) => {
    const s = v.section('chr17', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('bands', { label: 'Bands' }).layer('bands', { type: 'bands', data: bands });
    s.track('names', { label: 'Band' }).layer('ranges', { type: 'ranges', data: bands });
    v.navigator('chr17', ideogramNavigator(s, bands));
  });
  return b.build();
}

/** Generated BED12 genes: a block an exon, the thick part the CDS, the strand on each. */
export function genes(): Node {
  const b = createSeqSpec({ title: 'BED12 genes' });
  const sizes = b.data.source('sizes', {
    format: 'chrom-sizes',
    value: `${GeneratedSystem}\t${GeneratedLength}\n` as never,
  });
  const bed = b.data.source('genes-bed', { format: 'bed', uri: GeneratedBedFile, fields: { label: '3' } });
  const chr = b.data.sequence('chrG', { alphabet: 'dna', source: sizes, name: GeneratedSystem });
  const genes = b.data.ranges('genes', { target: chr, source: bed });
  b.views.view('main', (v) => {
    const s = v.section('genes', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('genes', { label: 'Genes' }).layer('ranges', { type: 'ranges', data: genes });
  });
  return b.build();
}

/** A misread BED: the lines read drawn, the rest said (a line number each) under the view. */
export function misread(): Node {
  const b = createSeqSpec({ title: 'A misread BED' });
  const sizes = b.data.source('sizes', { format: 'chrom-sizes', value: `${GeneratedSystem}\t7000\n` as never });
  const bed = b.data.source('misread-bed', { format: 'bed', uri: MisreadBedFile, fields: { label: '3' } });
  const chr = b.data.sequence('chrG', { alphabet: 'dna', source: sizes, name: GeneratedSystem });
  const rows = b.data.ranges('rows', { target: chr, source: bed });
  b.views.view('main', (v) => {
    const s = v.section('misread', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('rows', { label: 'Rows read' }).layer('ranges', { type: 'ranges', data: rows });
  });
  return b.build();
}

/**
 * What the views need besides the default registry: the `bed` runtime, its files read from their recordings where they
 * have one, each answer after `delayMs` (150) on `clock` (the synthetic twin the stories and tests share).
 */
export function registry(options: { readonly delayMs?: number; readonly clock?: StoreClock } = {}): SeqRegistryEntry[] {
  const runtime = delayedRuntime(withRecordedPages(BedRuntime), {
    name: 'bed-async',
    formats: BedFormats,
    delayMs: options.delayMs ?? 150,
    ...(options.clock ? { clock: options.clock } : {}),
  });
  return [{ kind: 'format-runtime', runtime }];
}
