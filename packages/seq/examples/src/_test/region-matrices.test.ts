/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type MatrixSummary, type MatrixItems, openSeq } from '@molstar/seqx-core/data';
import { describe, expect, it } from 'vitest';
import {
  GenotypeBin,
  GenotypeDeletion,
  GenotypeLength,
  GenotypeSamples,
  GenotypeSystem,
  genotypeItems,
  genotypesSource,
} from '../regions/genotypes.js';
import {
  MethylationBin,
  MethylationLength,
  MethylationSamples,
  MethylationSystem,
  methylationItems,
  methylationSource,
} from '../regions/methylation.js';

const signal = new AbortController().signal;

describe('the generated genotypes (Features/Region reads; EC-220, EC-222)', () => {
  it('are 40 samples over sorted variants, one long deletion across the SNVs inside it, some calls missing', () => {
    const g = genotypeItems();
    expect(g.states).toBe(40);
    expect(GenotypeSamples.length).toBe(40);
    for (let k = 1; k < g.start.length; k++) expect(g.start[k]).toBeGreaterThanOrEqual(g.start[k - 1]);
    expect([g.start[g.deletion], g.end[g.deletion]]).toEqual([GenotypeDeletion.start, GenotypeDeletion.end]);
    // the carriers have the deletion, the others none; SNVs inside it are the deletion’s to cover
    for (let s = 0; s < 40; s++)
      expect(g.codes![g.deletion * 40 + s]).toBe((GenotypeDeletion.carriers as readonly number[]).includes(s) ? 1 : 0);
    let inside = 0;
    for (let k = 0; k < g.start.length; k++)
      if (k !== g.deletion && g.start[k] > GenotypeDeletion.start && g.end[k] < GenotypeDeletion.end) inside++;
    expect(inside).toBeGreaterThan(20);
    expect(g.codes!.some((c) => c === 255)).toBe(true);
    // deterministic
    expect(Array.from(genotypeItems().codes!.subarray(0, 200))).toEqual(Array.from(g.codes!.subarray(0, 200)));
  });

  it('serve bins whose class is the items’ majority and whose share says how many agree', async () => {
    const src = openSeq(genotypesSource());
    const read = async (level: number) => {
      const a = await src.read(
        { system: GenotypeSystem, start: 0, end: GenotypeLength, level, shapes: ['matrix'] },
        signal,
      );
      if (!a.ok) throw new Error(a.issue.message);
      return a.tables.matrix!;
    };
    const items = (await read(1)) as MatrixItems;
    const bins = (await read(GenotypeBin)) as MatrixSummary;
    expect(bins.states).toBe(40);
    // a bin of one item agrees with it entirely
    const S = 40;
    const b = Math.floor(items.start[300] / GenotypeBin);
    let n = 0;
    for (let k = 0; k < items.start.length; k++) if (Math.floor(items.start[k] / GenotypeBin) === b) n++;
    if (n === 1) {
      for (let s = 0; s < S; s++) {
        const c = items.codes![300 * S + s];
        if (c !== 255) {
          expect(bins.classes![b * S + s]).toBe(c);
          expect(bins.shares![b * S + s]).toBe(255);
        }
      }
    }
    expect(bins.classes!.length).toBe(Math.ceil(GenotypeLength / GenotypeBin) * S);
  });
});

describe('the generated methylation (EC-221)', () => {
  it('has fractions in 0 to 1, a tenth missing as NaN, and an island of unmethylated CpGs', () => {
    const m = methylationItems();
    const v = m.value!;
    const nan = v.filter((x) => Number.isNaN(x)).length;
    expect(m.missing).toBe(nan);
    expect(nan / v.length).toBeGreaterThan(0.05);
    expect(nan / v.length).toBeLessThan(0.15);
    for (const x of v) if (!Number.isNaN(x)) expect(x >= 0 && x <= 1).toBe(true);
    const mean = (inIsland: boolean) => {
      let s = 0;
      let c = 0;
      for (let k = 0; k < m.start.length; k++) {
        const island = m.start[k] >= 240_000 && m.start[k] < 262_000;
        if (island !== inIsland) continue;
        for (let j = 0; j < m.states; j++) if (!Number.isNaN(v[k * m.states + j])) (s += v[k * m.states + j]), c++;
      }
      return s / c;
    };
    expect(mean(true)).toBeLessThan(0.15);
    expect(mean(false)).toBeGreaterThan(0.5);
  });

  it('serves bins of the mean with its counts, NaN where no CpG was covered', async () => {
    const src = openSeq(methylationSource());
    const a = await src.read(
      { system: MethylationSystem, start: 0, end: MethylationLength, level: MethylationBin, shapes: ['matrix'] },
      signal,
    );
    if (!a.ok) throw new Error(a.issue.message);
    const bins = a.tables.matrix as MatrixSummary;
    expect(bins.states).toBe(MethylationSamples.length);
    expect(bins.mean!.some((x) => Number.isNaN(x))).toBe(true);
    expect(bins.mean!.some((x) => !Number.isNaN(x))).toBe(true);
    for (let i = 0; i < bins.mean!.length; i++) expect(Number.isNaN(bins.mean![i])).toBe(bins.count![i] === 0);
  });
});
