/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { itWithin } from '@molstar/test-harness/timing';
import { expect } from 'vitest';
import { assembleTranscripts, parseGff } from '../gff.js';
import { transcriptTest } from '../gff-assemble.js';
import { GffCauses } from '../gff-issues.js';
import { generateGff3 } from './gff-fixtures.js';

// The GFF runtime's time (_spec/README.md "Data", GFF3 and GTF; performance.md's rows): 100,000 GFF3 lines assembled
// into transcripts under 500 ms. Its quick representative (10,000 lines, no clock) is in `gff.test.ts`.

itWithin(
  '100,000 GFF3 lines assembled under 500 ms',
  (within) => {
    const g = generateGff3(7000, ['chr1', 'chr2', 'chr3'], 4, 200_000_000);
    expect(g.lines).toBeGreaterThan(100_000);
    const file = parseGff(g.text, 'gff3', transcriptTest(), new GffCauses(), 1e7);
    if (typeof file === 'string') throw new Error(file);
    const a = within(500, () => assembleTranscripts(file.rows));
    expect(a.transcripts.length).toBe(g.transcripts);
  },
  60_000,
);
