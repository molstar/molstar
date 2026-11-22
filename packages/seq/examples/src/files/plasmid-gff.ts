/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Node, printAddress } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { PlasmidFiles } from '../plasmid/derive.js';
import { PlasmidGffFile, PlasmidSeqId } from './data.js';

// The plasmid's features from a GFF3 file, as a ring with a linear detail following it: the file's `region` row says the
// construct is circular; bla and a primer site end past its length (one row each, as NCBI writes a ring's features) and
// are drawn as one item across the origin. A track a feature type, a leading item form on `category`.

export const title = 'Plasmid features from GFF3';

/** The feature types of the file, a ring track each: the SO term and its label. */
export const GffTracks: readonly (readonly [type: string, label: string])[] = [
  ['CDS', 'CDS'],
  ['restriction_enzyme_recognition_site', 'Sites'],
  ['misc_feature', 'Misc'],
  ['rep_origin', 'Origin'],
  ['primer_bind', 'Primer sites'],
];

/** The ring and its detail; `uri` the file's location. */
export function build(uri: string = PlasmidGffFile): Node {
  const b = createSeqSpec({ title });
  // the letters from the rotated variant's `seq` file: once it is open its length is known, and a row past it wraps
  const record = b.data.source('record', { format: 'seq', uri: PlasmidFiles.rotated });
  const seq = b.data.sequence('plasmid', { alphabet: 'dna', source: record, cyclic: true, name: PlasmidSeqId });
  // GFF3 is a preset: nine columns, from 1, the last position named; the type is the category, `Name` the label
  const gff = b.data.source('features-file', { format: 'gff3', uri, fields: { category: 'type', label: 'Name' } });
  const features = b.data.ranges('features', { target: seq, source: gff });
  const data = (type: string) => `${features}${printAddress([{ kind: 'item', params: { category: type } }])}`;
  b.views.view('main', (v) => {
    const ring = v.section('ring', { of: seq, layout: 'rings' });
    for (const [type, label] of GffTracks)
      ring.track(type, { label }).layer('features', { type: 'ranges', data: data(type) });
    ring.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
    const detail = v.section('detail', { of: seq, follow: ring });
    for (const [type, label] of GffTracks.slice(0, 2))
      detail.track(type, { label }).layer('features', { type: 'ranges', data: data(type) });
    detail.track('letters', { label: PlasmidSeqId }).layer('letters', { type: 'letters' });
    detail.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}
