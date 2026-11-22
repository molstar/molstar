/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { itWithin } from '@molstar/test-harness/timing';
import { expect } from 'vitest';
import { DefaultOpenLimits } from '../../data/protocol.js';
import { VcfRuntime } from '../vcf.js';
import { memoryAssets } from './bed-fixtures.js';
import { generateVcf, indexedVcf } from './vcf-fixtures.js';

// The VCF runtime's time (_spec/README.md "Data", VCF; performance.md's row): a region of 10,000 records decoded under
// 100 ms, sites (its BGZF blocks inflated, its lines parsed, its rows answered). Its quick representative (1,000
// records) is in `vcf.test.ts`.

itWithin(
  'a VCF region of 10,000 records decoded under 100 ms (sites)',
  async (within) => {
    const g = generateVcf(10_000, ['chr1'], 13, 1_000_000);
    const { file, index } = indexedVcf(g, 65_280);
    const signal = new AbortController().signal;
    const s = await VcfRuntime.open(
      { address: '/data/v', format: 'vcf', uri: 'f', index: 'i' },
      { assets: memoryAssets({ f: file, i: index }), signal, limits: DefaultOpenLimits },
    );
    const a = await within(100, () =>
      s.read({ system: 'chr1', start: 0, end: 1_100_000, level: 1, shapes: ['ranges'] }, signal),
    );
    expect(a.ok && a.tables.ranges?.form === 'list' && a.tables.ranges.start.length).toBeGreaterThanOrEqual(10_000);
    s.dispose();
  },
  60_000,
);
