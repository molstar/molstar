/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';

// A genome: an assembly from its sequences, and a genome section of it (a block a sequence, a ruler in bases); the
// second section's `include` shows one chromosome (a locus is a query, not a builder call). A link over the assembly
// names its ends' sequences, `a_name` and `b_name`, as a map's runs name theirs.

export function build(): Node {
  const b = createSeqSpec({ title: 'A genome view' });
  const toy = b.data.assembly(
    'toy',
    assembly('Toy1', [
      { name: 'chr1', length: 48_000_000 },
      { name: 'chr2', length: 36_000_000 },
      { name: 'chr3', length: 22_000_000 },
      { name: 'chrM', length: 16_569, cyclic: true },
    ]),
  );
  const link = b.data.pairs('link', {
    a: toy,
    rows: {
      a_start: [12_000_000],
      a_end: [14_000_000],
      a_name: ['chr1'],
      b_start: [5_000_000],
      b_end: [7_000_000],
      b_name: ['chr3'],
    },
  });
  b.views.view('main', (v) => {
    genome(v, 'all', toy, {}, (s) => s.track('link').layer('link', { type: 'pairs', data: link, legend: 'none' }));
    genome(v, 'locus', toy, { include: ['chr2'] });
  });
  return b.build();
}
