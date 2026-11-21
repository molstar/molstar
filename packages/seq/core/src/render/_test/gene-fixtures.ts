/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { generateGff3 } from '../../runtimes/_test/gff-fixtures.js';

// Gene models for the frame scenarios (`scenarios.ts`): a generated GFF3 of 400 genes over 2 Mb read by the GFF runtime,
// its transcripts drawn as the ranges layer's gene line, the sequence's length from its `##sequence-region`.

let text: string | undefined;

export function geneModelsSpec(): Node {
  text ??= generateGff3(400, ['chr1'], 7).text;
  const b = createSeqSpec({ title: 'gene models' });
  const file = b.data.source('gff', { format: 'gff3', value: text as never });
  const chr = b.data.sequence('chr1', { alphabet: 'dna', source: file, name: 'chr1' });
  const tx = b.data.ranges('tx', { target: chr, source: file });
  b.views.view('main', (v) => {
    const s = v.section('chr1', { of: chr });
    s.track('ruler', { label: '' }).layer('ruler', { type: 'ruler', units: 'bases' });
    s.track('tx', { label: 'Transcripts' }).layer('tx', { type: 'ranges', data: tx, line: 'gene' } as never);
  });
  return b.build();
}
