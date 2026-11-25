/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The places index's allocation (_spec/README.md "Resolved view", Places): a query into a caller's arrays (what a
// layer reads per arc or ribbon end) allocates nothing once warm. By hand: `pnpm -C packages/seq/core test:perf`.

import { createSeqSpec } from '@molstar/seqx-core/lang';
import { resolveSpec } from '@molstar/seqx-core/resolve';
import { quickCase } from '@molstar/test-harness/quick';
import { describe, expect, inject, it } from 'vitest';
import { resolveView } from '../index.js';
import { placesInto } from '../places.js';

describe('the places index', () => {
  it('a query over 3,000 blocks into the caller’s arrays: 0 bytes warm', ({ task }) => {
    const b = createSeqSpec();
    const s = b.data.sequence('s', { alphabet: 'dna', letters: 'A'.repeat(3000 * 700 + 1000) });
    b.views.view('v', (v) =>
      v.section('a', {}, (sec) => {
        for (let i = 0; i < 3000; i++) sec.block(`b${i}`, { of: s, range: [i * 700, i * 700 + 1000] } as never);
      }),
    );
    const section = resolveView(resolveSpec(b.build()), '/views/v')!.sections[0];
    const blocks = new Int32Array(8);
    const units = new Float64Array(8);
    const { r } = quickCase(task, 'places query', (i) => {
      // a whole position (an arc end): a fractional number passed to a call is the caller's box
      placesInto(section, '/data/s', (i * 7013) % 2_100_000, blocks, units);
    });
    console.info(`places query: ${r.bytesPerStep.toFixed(0)} bytes a query (${r.note})`);
    expect(r.bytesPerStep, r.note).toBeLessThanOrEqual(inject('perfThreshold'));
  });
});
