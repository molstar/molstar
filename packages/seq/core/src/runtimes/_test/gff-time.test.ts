/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// @module-tag heavy

import { itWithin } from '@molstar/test-harness/timing';
import { expect } from 'vitest';
import { DefaultOpenLimits } from '../../data/protocol.js';
import { assembleTranscripts, GffRuntime, parseGff } from '../gff.js';
import { transcriptTest } from '../gff-assemble.js';
import { GffCauses } from '../gff-issues.js';
import { memoryAssets } from './bed-fixtures.js';
import { generateGff3 } from './gff-fixtures.js';

// The GFF runtime's time (_spec/README.md "Data", GFF3 and GTF; performance.md's rows): 100,000 GFF3 lines assembled
// into transcripts under 500 ms (its quick representative, 10,000 lines and no clock, in `gff.test.ts`); `features:
// rows` of 40,000 lines of one `ID` opened under 500 ms (its pieces collected once: 65 s when each row merged them
// again); 2,000 region reads beside a chromosome-long feature under 200 ms (the long feature apart, not scanned from);
// 50 lines of a million tabs each refused under 200 ms (the tabs counted, the line never split whole).

const signal = new AbortController().signal;
const ctx = { assets: memoryAssets({}), signal, limits: DefaultOpenLimits };

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

itWithin(
  '`features: rows`: 40,000 lines of one ID opened under 500 ms',
  async (within) => {
    const lines: string[] = [];
    for (let e = 0; e < 40_000; e++) lines.push(`c\t.\tCDS\t${e * 10 + 1}\t${e * 10 + 5}\t.\t+\t0\tID=cds1`);
    const text = lines.join('\n');
    const s = await within(500, () =>
      GffRuntime.open({ address: '/s', format: 'gff3', value: text as never, features: 'rows' } as never, ctx),
    );
    expect(s.coverage.systems.length).toBe(1);
  },
  120_000,
);

itWithin(
  '2,000 region reads beside a chromosome-long feature under 200 ms',
  async (within) => {
    const lines = ['c\t.\trepeat_region\t1\t100000000\t.\t+\t.\tID=long'];
    for (let i = 0; i < 100_000; i++)
      lines.push(`c\t.\trepeat_region\t${i * 1000 + 1}\t${i * 1000 + 100}\t.\t+\t.\tID=r${i}`);
    const s = await GffRuntime.open({ address: '/s', format: 'gff3', value: lines.join('\n') as never }, ctx);
    const read = async () => {
      for (let k = 0; k < 2000; k++) {
        const start = 50_000_000 + k * 1000;
        const a = await s.read({ system: 'c', start, end: start + 1000, level: 1, shapes: ['ranges'] }, signal);
        if (!a.ok) throw new Error('not read');
      }
    };
    await within(200, read);
  },
  120_000,
);

itWithin(
  '50 lines of a million tabs each refused under 200 ms',
  (within) => {
    const line = `c\t.\tgene\t1\t10\t.\t+\t.\tID=a${'\t'.repeat(1_000_000)}`;
    const text = Array.from({ length: 50 }, () => line).join('\n');
    const causes = new GffCauses();
    const f = within(200, () => parseGff(text, 'gff3', transcriptTest(), causes, 1e6));
    expect(typeof f === 'string' ? f : f.rows.length).toBe(0);
  },
  120_000,
);
