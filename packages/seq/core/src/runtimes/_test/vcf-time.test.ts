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
import { flagsLine, generateVcf, indexedVcf, tablessLines, VcfHeaderText } from './vcf-fixtures.js';

// The VCF runtime's time (_spec/README.md "Data", VCF): a region of 10,000 records decoded under
// 100 ms, sites (its BGZF blocks inflated, its lines parsed, its rows answered), and with 100 samples' genotypes under
// 200 ms. Their quick representatives (1,000 records) are in `vcf.test.ts`. A hostile file's lines cost their own length: 40,000 lines without tabs, and a line of
// 100,000 flags, read whole at once (their quick counts in `vcf-cost.test.ts`).

/** `body` after the header, read whole by the runtime. */
async function openWhole(body: Uint8Array) {
  const head = new TextEncoder().encode(VcfHeaderText);
  const text = new Uint8Array(head.length + body.length);
  text.set(head);
  text.set(body, head.length);
  return VcfRuntime.open(
    { address: '/data/v', format: 'vcf', uri: 'f' },
    { assets: memoryAssets({ f: text }), signal: new AbortController().signal, limits: DefaultOpenLimits },
  );
}

itWithin('40,000 lines without their tabs read whole under 500 ms (each search in its line)', async (within) => {
  const body = tablessLines(40_000);
  const s = await within(500, () => openWhole(body));
  expect(s.issues?.[0].message).toMatch(/^Fewer than eight columns: 40000/);
});

itWithin('a line of 100,000 INFO flags read whole under 200 ms (`=` looked for in each field)', async (within) => {
  const body = flagsLine(100_000);
  const s = await within(200, () => openWhole(body));
  expect(s.coverage.systems.map((x) => x.name)).toEqual(['chr1']);
});

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

itWithin(
  'a VCF region of 10,000 records with 100 samples’ genotypes decoded under 200 ms (the matrix’s codes)',
  async (within) => {
    const g = generateVcf(10_000, ['chr1'], 14, 1_000_000, { samples: 100, genotypes: true });
    const { file, index } = indexedVcf(g, 65_280);
    const signal = new AbortController().signal;
    const s = await VcfRuntime.open(
      { address: '/data/v', format: 'vcf', uri: 'f', index: 'i' },
      { assets: memoryAssets({ f: file, i: index }), signal, limits: DefaultOpenLimits },
    );
    const a = await within(200, () =>
      s.read({ system: 'chr1', start: 0, end: 1_100_000, level: 1, shapes: ['ranges', 'matrix'] }, signal),
    );
    expect(a.ok && a.tables.matrix?.form === 'items' && a.tables.matrix.codes!.length).toBeGreaterThanOrEqual(
      10_000 * 100,
    );
    s.dispose();
  },
  60_000,
);
