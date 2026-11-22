/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';

// A genome: an assembly from its sequences, and a genome section of it (a block a sequence, a ruler in bases); the
// second view's `include` is its locus, the chromosomes it shows.

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
  b.views.view('main', (v) => {
    genome(v, 'all', toy);
    genome(v, 'locus', toy, { include: ['chr2'] });
  });
  return b.build();
}
