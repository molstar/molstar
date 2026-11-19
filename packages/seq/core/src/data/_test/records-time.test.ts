/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

// The open's target (_spec/README.md "Data", Records): 100,000 TSV rows read into the protocol's tables under 300 ms,
// the best of three (`itWithin`). The quick lane holds a representative of 10,000 rows (`records.test.ts`).

import { itWithin } from '@molstar/test-harness/timing';
import { describe, expect } from 'vitest';
import { DefaultOpenLimits, type RangesList, type RegionAnswer } from '../protocol.js';
import { readRecords } from '../records.js';

function table(n: number): string {
  const lines = ['chrom\tstart\tend\tstrand\tname\tscore'];
  for (let i = 0; i < n; i++) {
    const s = ((i * 7919) % 1_000_000) + 1;
    lines.push(`chr${i % 3}\t${s}\t${s + 50}\t${i % 2 ? '+' : '-'}\tf${i}\t${i % 100}`);
  }
  return `${lines.join('\n')}\n`;
}

describe('records: the open', () => {
  itWithin('reads 100,000 TSV rows under 300 ms', (within) => {
    const text = table(100_000);
    const r = within(300, () =>
      readRecords(
        { format: 'tsv', text },
        { fields: { target: 'chrom' }, types: { score: 'float32' } },
        { system: 's', limits: DefaultOpenLimits },
      ),
    );
    let rows = 0;
    for (const p of r.data!.pairs)
      rows += ((p.answer as Extract<RegionAnswer, { ok: true }>).tables.ranges as RangesList).start.length;
    expect(rows).toBe(100_000);
  });
});
