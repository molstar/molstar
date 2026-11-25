/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type GridDense, openSeq, type SeqData } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import { HicBin, HicChromosomes, HicCoarseBin, hicOneLevel, hicSource } from '../regions/hic.js';

const signal = new AbortController().signal;
const total = (v: Float32Array) => v.reduce((a, x) => a + x, 0);

async function tile(data: SeqData, x: string, y: string, level: number): Promise<GridDense> {
  const src = openSeq(data);
  const cx = HicChromosomes.find((c) => c.name === x)!;
  const cy = HicChromosomes.find((c) => c.name === y)!;
  const a = await src.read(
    { system: x, start: 0, end: cx.length, level, shapes: ['grid'], y: { system: y, start: 0, end: cy.length } },
    signal,
  );
  if (!a.ok) throw new Error(a.issue.message);
  return a.tables.grid as GridDense;
}

describe('the generated Hi-C (Features/Region reads)', () => {
  it('is deterministic, its upper triangle of block pairs served, symmetric within a chromosome', async () => {
    const a = hicSource();
    const b = hicSource();
    const x = await tile(a, 'chr5', 'chr5', HicBin);
    expect(total((await tile(b, 'chr5', 'chr5', HicBin)).value)).toBe(total(x.value));
    // a chromosome with itself is symmetric; the diagonal the strongest
    for (let i = 0; i < x.nx; i += 7)
      for (let j = 0; j < x.ny; j += 11) expect(x.value[i * x.ny + j]).toBe(x.value[j * x.ny + i]);
    expect(x.value[3 * x.ny + 3]).toBeGreaterThan(x.value[3 * x.ny + 60]);
    // the blocks below the diagonal are not served (the reader transposes)
    expect(a.coverage.shapes.find((s) => s.system === 'chr2')!.y).toEqual(['chr2', 'chr3', 'chr4', 'chr5', 'chrX']);
    expect(a.coverage.shapes.find((s) => s.system === 'chrX')!.y).toEqual(['chrX']);
  });

  it('has two levels summing to the finest, a whole-assembly system listed, and a one-level twin of the same counts', async () => {
    const two = hicSource();
    expect(two.coverage.levels).toEqual([HicBin, HicCoarseBin]);
    expect(hicOneLevel().coverage.levels).toEqual([HicBin]);
    const length = HicChromosomes.reduce((a, c) => a + c.length, 0);
    expect(two.coverage.systems.at(-1)).toEqual({ name: 'all', length });
    expect(two.coverage.shapes[0].whole).toEqual({ system: 'all', sequences: HicChromosomes.map((c) => c.name) });
    const fine = await tile(two, 'chr3', 'chr4', HicBin);
    const coarse = await tile(two, 'chr3', 'chr4', HicCoarseBin);
    expect(total(coarse.value)).toBe(total(fine.value));
    expect(coarse.nx).toBe(Math.ceil(fine.nx / 4));
    expect(total((await tile(hicOneLevel(), 'chr3', 'chr4', HicBin)).value)).toBe(total(fine.value));
  });

  it('has a translocation’s hot spot between chr2 and chr5 and noise elsewhere', async () => {
    const two = hicSource();
    const hot = await tile(two, 'chr2', 'chr5', HicBin);
    const cold = await tile(two, 'chr3', 'chr4', HicBin);
    expect(Math.max(...hot.value)).toBeGreaterThan(20);
    expect(Math.max(...cold.value)).toBeLessThanOrEqual(2);
  });
});
