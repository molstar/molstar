/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The ring window's fixture (CG3): a cyclic sequence on rings and a linear section following it, for the tests and the
// heap guards.

import { createSeqViewSpec } from '@molstar/seq-core/builder';

export const N = 600;

/** What `ringWindowDoc` varies: more linear followers, the detail's piece, a second (linear) block before the ring's. */
export interface RingWindowOptions {
  readonly followers?: number;
  readonly range?: readonly [number, number];
  readonly before?: boolean;
}

/**
 * A cyclic sequence of `n` nucleotides (600) on a ring (ruler, features) and a linear section `lin` following it (with
 * `followers`, `lin2`… too); the ring's origin `origin`.
 */
export function ringWindowDoc(origin = 0, n = N, o: RingWindowOptions = {}) {
  const b = createSeqViewSpec();
  const s = b.data.sequence({ ref: 's', value: 'ACGTACGTTA'.repeat(n / 10), alphabet: 'dna', cyclic: true });
  const f = b.data.annotation({
    ref: 'f',
    geometry: 'ranges',
    target: s,
    value: {
      start: [n / 6, n / 2, (11 * n) / 12],
      stop: [n / 3, (2 * n) / 3, n / 12],
      label: ['alpha', 'beta', 'gamma'],
    },
  } as never);
  const v = b.view({ ref: 'v', label: 'V' });
  const ring = v.section({ ref: 'ring', style: 'circular' });
  if (o.before) {
    const s2 = b.data.sequence({ ref: 's2', value: 'ACGTACGTTA'.repeat(n / 10), alphabet: 'dna' });
    ring.block({ ref: 'ring.b2', sequence: s2 });
  }
  ring.block({ ref: 'ring.b', sequence: s, origin });
  ring.track({ label: '' }).ruler();
  ring.track({ label: 'F', height: 2 }).ranges({ data: f, labels: true } as never);
  for (let k = 1; k <= (o.followers ?? 1); k++) {
    const lin = v
      .section({ ref: k === 1 ? 'lin' : `lin${k}`, follow: ring.handle })
      .block({ sequence: s, ...(o.range ? { range: [...o.range] } : {}) });
    lin.track({ label: '' }).ruler();
    lin.track({ label: 'F' }).ranges({ data: f } as never);
  }
  return b.getDocument();
}
