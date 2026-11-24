/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { P53 } from '../p53/data.js';
import { Cdk2, type DomainFile, DomainFiles, FilesFolder, ManySequencesFile, ManySystems } from './data.js';

// p53's domains read from a table file: the same six rows as TSV, CSV, JSON and MessagePack, an array of records or an
// object of columns, each drawn the same. The tables count from 1 and name their last position (the file defaults), so
// the source says nothing but its format and its location.

/** The view of one file: the protein's letters and the domains the file holds. */
export function build(file: DomainFile = 'domains.tsv'): Node {
  const b = createSeqSpec({ title: `p53 domains from ${file}` });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const domains = b.data.source('domains-file', {
    format: DomainFiles[file].format,
    uri: `${FilesFolder}${file}`,
    // the file's rows are about this sequence: the system the source's own rows are on
    system: 'p53',
  });
  const ranges = b.data.ranges('domains', { target: p53, source: domains });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    s.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: ranges });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
  });
  return b.build();
}

// One file, two sequences: the rows name theirs by a `target` column, one row has none (the source's own system, p53).

/** The view of the many-sequences table: a section a sequence, each showing only its own rows. */
export function buildMany(file: string = `${FilesFolder}${ManySequencesFile}`): Node {
  const b = createSeqSpec({ title: 'Two sequences from one table' });
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const cdk2 = b.data.sequence('cdk2', { alphabet: 'protein', letters: Cdk2 });
  // `target` names each row's sequence (a sequence reads the system of its name); the empty one is on `system`
  const rows = b.data.source('many-file', { format: 'tsv', uri: file, system: ManySystems.own });
  const p53Ranges = b.data.ranges('p53-features', { target: p53, source: rows });
  const cdk2Ranges = b.data.ranges('cdk2-features', { target: cdk2, source: rows });
  b.views.view('main', (v) => {
    for (const [id, label, of, ranges] of [
      ['p53', 'p53', p53, p53Ranges],
      ['cdk2', 'CDK2', cdk2, cdk2Ranges],
    ] as const) {
      const s = v.section(id, { of });
      s.track('letters', { label }).layer('letters', { type: 'letters' });
      s.track('features', { label: 'Features' }).layer('ranges', { type: 'ranges', data: ranges });
      s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler' });
    }
  });
  return b.build();
}
