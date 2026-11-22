/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, printAddress } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import { type AssetReader, DefaultOpenLimits, delayedRuntime, type StoreClock } from '@molstar/seqx-core/data';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { GenbankFormats, GenbankRuntime, type GenbankSource } from '@molstar/seqx-core/runtimes/genbank';
import { GffTracks } from '../files/plasmid-gff.js';
import { PlasmidFiles, PlasmidName } from '../plasmid/derive.js';
import { MisreadGenbankFile, MisreadGenbankSystem } from './files.js';

// GenBank's views (Features/Genome formats/GenBank): the plasmid's generated record (`fixtures/plasmid/
// pUC19c-rotated.gb`) read by the opt-in `genbank` runtime, its letters and its features from the one file, drawn as
// the plasmid's GFF3 is (Features/Files: the same ring and detail, the same tracks by feature key); and a misread file
// and its issues. The record says it is circular; bla's `complement(join(2312..2686,1..486))` is one item across the
// origin, as the GFF3's row ending past the length is.

/** The sections' addresses (what the tests read). */
export const PlasmidRingSection = '/views/main/ring';
export const MisreadGenbankSection = '/views/main/misread';
/** The features' annotation (an item's address is `<it>/(item :id "CDS@2312-3172")`). */
export const FeaturesAddress = '/data/features';

/** The plasmid from its GenBank record: the GFF3 view's ring and detail, letters and features from the one file. */
export function plasmid(uri: string = PlasmidFiles.rotatedGenbank): Node {
  const b = createSeqSpec({ title: 'Plasmid features from GenBank' });
  const gb = b.data.source('record', { format: 'genbank', uri });
  const seq = b.data.sequence('plasmid', { alphabet: 'dna', source: gb, cyclic: true, name: PlasmidName });
  const features = b.data.ranges('features', { target: seq, source: gb });
  // GenBank's key for the GFF3's restriction sites is `protein_bind`
  const data = (key: string) =>
    `${features}${printAddress([{ kind: 'item', params: { category: key === 'restriction_enzyme_recognition_site' ? 'protein_bind' : key } }])}`;
  b.views.view('main', (v) => {
    const ring = v.section('ring', { of: seq, layout: 'rings' });
    for (const [key, label] of GffTracks)
      ring.track(key, { label }).layer('features', { type: 'ranges', data: data(key) });
    ring.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
    const detail = v.section('detail', { of: seq, follow: ring });
    for (const [key, label] of GffTracks.slice(0, 2))
      detail.track(key, { label }).layer('features', { type: 'ranges', data: data(key) });
    detail.track('letters', { label: PlasmidName }).layer('letters', { type: 'letters' });
    detail.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}

/** A misread GenBank file: what reads drawn, the rest said (a line number each) under the view. */
export function misread(): Node {
  const b = createSeqSpec({ title: 'A misread GenBank file' });
  const file = b.data.source('misread-gb', { format: 'genbank', uri: MisreadGenbankFile });
  const seq = b.data.sequence('misread', { alphabet: 'dna', source: file, name: MisreadGenbankSystem });
  const features = b.data.ranges('features', { target: seq, source: file });
  b.views.view('main', (v) => {
    const s = v.section('misread', { of: seq });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('features', { label: 'Read' }).layer('features', { type: 'ranges', data: features });
  });
  return b.build();
}

/**
 * A GenBank file's records as a few lines (name, accession, length, topology, molecule, DEFINITION, organism,
 * features), read as a view reads the file, by the `genbank` runtime: what a story prints under the view.
 */
export async function genbankSummary(uri: string, assets: AssetReader, signal: AbortSignal): Promise<string> {
  const s = (await GenbankRuntime.open(
    { address: '/data/file', format: 'genbank', uri },
    { assets, signal, limits: DefaultOpenLimits },
  )) as GenbankSource;
  try {
    return s.records
      .map((r) =>
        [
          `${r.name}${r.accession ? ` (${r.accession})` : ''}: ${r.length} ${r.molecule === 'AA' ? 'aa' : 'bp'}, ${r.circular ? 'circular' : 'linear'}, ${r.molecule}`,
          `  ${r.definition ?? '(no DEFINITION)'}`,
          `  ${r.organism ?? '(no ORGANISM)'}; ${r.features} features`,
        ].join('\n'),
      )
      .join('\n');
  } finally {
    s.dispose();
  }
}

/** What the views need besides the default registry: the `genbank` runtime, each answer after `delayMs` (150) on `clock`. */
export function registry(options: { readonly delayMs?: number; readonly clock?: StoreClock } = {}): SeqRegistryEntry[] {
  const timing = { delayMs: options.delayMs ?? 150, ...(options.clock ? { clock: options.clock } : {}) };
  return [
    {
      kind: 'format-runtime',
      runtime: delayedRuntime(GenbankRuntime, { name: 'genbank-async', formats: GenbankFormats, ...timing }),
    },
  ];
}
