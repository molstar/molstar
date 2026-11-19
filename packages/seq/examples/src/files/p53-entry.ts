/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { uniprotEntry } from './uniprot.js';

// p53 from UniProt's own file: the REST entry of P04637 (`fixtures/p53/P04637.json`, trimmed, never reshaped), its
// letters and its features read by the records runtime, a track a feature type (one read of the file, a layer a type
// with a leading item form).

export const title = 'p53 from UniProt’s entry';

/** The entry's location: the fixture, or the service's URL in a showcase. */
export const EntryFile = 'fixtures/p53/P04637.json';

/** The feature types the view shows, a track each, with the track's label. */
export const EntryTracks: readonly (readonly [type: string, label: string])[] = [
  ['DNA binding', 'DNA binding'],
  ['Region', 'Regions'],
  ['Motif', 'Motifs'],
  ['Helix', 'Helices'],
  ['Beta strand', 'Strands'],
  ['Turn', 'Turns'],
  ['Modified residue', 'Modified residues'],
];

/** The view's spec; `uri` the entry's location. */
export function build(uri: string = EntryFile): Node {
  const b = createSeqSpec({ title });
  const entry = uniprotEntry(b, uri, { accession: 'P04637' });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: entry.protein });
    s.track('letters', { label: 'P04637' }).layer('letters', { type: 'letters' });
    for (const [type, label] of EntryTracks)
      s.track(type.toLowerCase().replaceAll(' ', '-'), { label }).layer('ranges', {
        type: 'ranges',
        data: entry.data(type),
      });
    s.track('ruler', { label: 'Ruler' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}
