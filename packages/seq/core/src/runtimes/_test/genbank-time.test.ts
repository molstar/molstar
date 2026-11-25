/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { itWithin } from '@molstar/test-harness/timing';
import { expect } from 'vitest';
import { DefaultOpenLimits } from '../../data/protocol.js';
import { GenbankRuntime } from '../genbank.js';
import { memoryAssets } from './bed-fixtures.js';
import { generated } from './genbank-fixtures.js';

// The GenBank runtime's time (_spec/README.md "Data", GenBank): 10,000 features opened (parsed, their locations read,
// held by sequence, the name index built) under 100 ms (its quick representative, 1,000 features and no clock, in
// `genbank.test.ts`); a file of 50 lines of a million spaces each and a quote never closed under 200 ms.

const signal = new AbortController().signal;
const ctx = { assets: memoryAssets({}), signal, limits: DefaultOpenLimits };

itWithin(
  '10,000 features parsed under 100 ms',
  async (within) => {
    const text = generated(10_000);
    const s = await within(100, () =>
      GenbankRuntime.open({ address: '/g', format: 'genbank', value: text as never }, ctx),
    );
    expect((s.issues ?? []).length).toBe(0);
  },
  60_000,
);

itWithin(
  'lines of a million spaces and a quote never closed under 200 ms',
  async (within) => {
    const wide = `                     /note="${' '.repeat(1_000_000)}`;
    const text = [
      'LOCUS       s   4 bp    DNA     linear',
      'FEATURES             Location/Qualifiers',
      '     misc_feature    1..2',
      ...Array.from({ length: 50 }, () => wide),
      'ORIGIN',
      '        1 acgt',
      '//',
    ].join('\n');
    const s = await within(200, () =>
      GenbankRuntime.open({ address: '/g', format: 'genbank', value: text as never }, ctx),
    );
    expect(s.coverage.systems[0].length).toBe(4);
  },
  60_000,
);
