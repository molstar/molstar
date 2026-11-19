/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { P53 } from '../p53/data.js';
import { type DomainFile, DomainFiles, FilesFolder } from './data.js';

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
