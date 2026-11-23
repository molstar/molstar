/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { itWithin } from '@molstar/test-harness/timing';
import { expect } from 'vitest';
import { DefaultOpenLimits } from '../../data/protocol.js';
import { BedRuntime } from '../bed.js';
import { indexPairs, PairsBuilder, pairsIn, pairsSchema } from '../bed-pairs.js';
import { BedCauses } from '../bed-issues.js';
import { generateBedpe, memoryAssets } from './bed-fixtures.js';

// The BED runtime's times (_spec/README.md "Data", BED family): 600,000 BEDPE rows, the most
// `maxPairsBytes` holds of a six-column BEDPE, indexed by both ends under 300 ms; a chrom.sizes of a million lines read
// in one pass. Their quick representatives (12,000 links, 20,000 lines) are in `bed.test.ts`.

const names = Array.from({ length: 24 }, (_, i) => `chr${i + 1}`);

itWithin(
  '600,000 BEDPE rows indexed by both ends under 300 ms',
  (within) => {
    const schema = pairsSchema(6);
    if ('error' in schema) throw new Error(schema.error);
    const b = new PairsBuilder(schema);
    const causes = new BedCauses();
    const text = generateBedpe(600_000, names, 1, 200_000_000, false);
    for (const l of text.trim().split('\n')) b.add(l.split('\t').slice(0, 6), '', causes);
    const rows = b.finish();
    const index = within(300, () => indexPairs(rows));
    expect(pairsIn(rows, index, 0, 0, 1_000_000).length).toBeGreaterThan(0);
  },
  60_000,
);

itWithin(
  'a chrom.sizes of 1,000,000 lines read under 1.5 s (a name seen kept in a set: linear)',
  async (within) => {
    const text = Array.from({ length: 1_000_000 }, (_, i) => `contig${i}\t${1000 + (i % 5000)}`).join('\n');
    const s = await within(1500, () =>
      BedRuntime.open(
        { address: '/data/sizes', format: 'chrom-sizes', value: text },
        { assets: memoryAssets({}), signal: new AbortController().signal, limits: DefaultOpenLimits },
      ),
    );
    expect(s.coverage.systems.length).toBe(1_000_000);
  },
  120_000,
);
