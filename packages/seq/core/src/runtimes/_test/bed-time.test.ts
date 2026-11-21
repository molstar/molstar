/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { itWithin } from '@molstar/test-harness/timing';
import { expect } from 'vitest';
import { indexPairs, PairsBuilder, pairsIn, pairsSchema } from '../bed-pairs.js';
import { BedCauses } from '../bed-issues.js';
import { generateBedpe } from './bed-fixtures.js';

// The BEDPE index's time (_spec/README.md "Data", BED family; performance.md's row): 600,000 rows, the most
// `maxPairsBytes` holds of a six-column BEDPE, indexed by both ends under 300 ms. The quick representative (12,000
// links, the milestone) is in `bed.test.ts`.

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
