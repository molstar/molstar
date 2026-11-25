/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { MatrixItems, SeqData, SeqPair } from '@molstar/seqx-core/data';
import { binValues } from '../../regions/methylation.js';
import type { ToyAssembly } from './assemblies.js';
import type { ToyGene } from './genes.js';
import { Concern, stream } from './streams.js';

// Step 6 of the Genomes generator (EC-220, EC-221): methylation of 20 samples at 250 CpGs of toy
// A's chr1 (A's own `CG`s in a 50 kb stretch, an island where they fall to near zero, 8 % of the values missing: NaN),
// and the same layer over the 150 genes (a gene-by-sample heatmap: each gene's mean over its body, a group of samples
// lower on a set of genes). Both in the `seq` format's matrix (the wide form): items at level 1, each bin's mean at
// 1,024 or 4,096 (what a view reads at the fit), written in MessagePack (float32 columns as bytes). The long form (a row a sample × position, pivoted by a records transform at
// open) is not written: no records transform pivots a table to a matrix yet.

export const MethylationSamples: readonly string[] = Array.from(
  { length: 20 },
  (_, i) => `${['liver', 'brain', 'blood', 'lung'][i % 4]}-${Math.floor(i / 4) + 1}`,
);
export const CpgCount = 250;
/** The stretch of chr1 the CpGs are in, and the island inside it. */
export const CpgStretch: readonly [number, number] = [220_000, 270_000];
export const CpgIsland: readonly [number, number] = [240_000, 246_000];
/** The share of values missing. */
export const MissingShare = 0.08;
/** The coarser level: each bin's mean (the CpGs' in 1,024 positions; the genes' in 4,096, a gene's slot). */
export const MethylationBin = 1024;
export const GeneMethylationBin = 4096;

/** `n` CpGs of `chr`'s letters in `[s, e)`, evenly picked among its `CG`s (a position the C's). */
function cpgs(letters: string, s: number, e: number, n: number): number[] {
  const all: number[] = [];
  for (let p = s; p < e - 1; p++) if (letters[p] === 'C' && letters[p + 1] === 'G') all.push(p);
  return Array.from({ length: Math.min(n, all.length) }, (_, i) => all[Math.floor((i * all.length) / n)]);
}

/** The matrix's two levels on each of `systems`' items (sequences without items left out). */
function levels(a: ToyAssembly, items: ReadonlyMap<string, MatrixItems>, bin: number): SeqData {
  const pairs: SeqPair[] = [];
  const systems = a.sequences.filter((s) => items.has(s.name));
  for (const s of systems) {
    const m = items.get(s.name)!;
    const request = (level: number) => ({
      system: s.name,
      start: 0,
      end: s.length,
      level,
      shapes: ['matrix' as const],
    });
    pairs.push({ request: request(1), answer: { ok: true, tables: { matrix: m } } });
    pairs.push({
      request: request(bin),
      answer: { ok: true, tables: { matrix: binValues(m, s.length, bin) } },
    });
  }
  return {
    coverage: {
      systems: systems.map((s) => ({ name: s.name, length: s.length })),
      shapes: systems.map((s) => ({ system: s.name, shape: 'matrix' as const, states: [...MethylationSamples] })),
      levels: [1, bin],
    },
    pairs,
  };
}

/** The CpGs' methylation (a fraction per CpG and sample, NaN missing), and how many are missing. */
export function cpgMethylation(
  a: ToyAssembly,
  letters: ReadonlyMap<string, string>,
  seed: number,
): { data: SeqData; positions: readonly number[]; missing: number } {
  const S = MethylationSamples.length;
  const pos = cpgs(letters.get('chr1')!, CpgStretch[0], CpgStretch[1], CpgCount);
  const value = new Float32Array(pos.length * S);
  let missing = 0;
  pos.forEach((p, i) => {
    const rng = stream(seed, Concern.methylation, i);
    const island = p >= CpgIsland[0] && p < CpgIsland[1];
    for (let s = 0; s < S; s++) {
      // brain samples lower outside the island: a tissue's mark
      const base = island ? 0.06 : s % 4 === 1 ? 0.55 : 0.82;
      const v = Math.min(1, Math.max(0, base + (rng.float() - 0.5) * 0.25));
      const gone = rng.float() < MissingShare;
      if (gone) missing++;
      value[i * S + s] = gone ? Number.NaN : Math.round(v * 1000) / 1000;
    }
  });
  const items: MatrixItems = {
    form: 'items',
    start: Int32Array.from(pos),
    end: Int32Array.from(pos, (p) => p + 2),
    states: S,
    value,
  };
  return { data: levels(a, new Map([['chr1', items]]), MethylationBin), positions: pos, missing };
}

/** The genes' methylation: each gene's mean over its body a sample (a quarter of the genes lower in blood). */
export function geneMethylation(a: ToyAssembly, genes: readonly ToyGene[], seed: number): SeqData {
  const S = MethylationSamples.length;
  const by = new Map<string, ToyGene[]>();
  for (const g of genes) by.set(g.chr, [...(by.get(g.chr) ?? []), g]);
  const items = new Map<string, MatrixItems>();
  for (const [chr, list] of by) {
    const value = new Float32Array(list.length * S);
    list.forEach((g, i) => {
      const rng = stream(seed, Concern.methylation, 10_000 + genes.indexOf(g));
      const level = 0.2 + 0.6 * rng.float();
      const marked = rng.float() < 0.25;
      for (let s = 0; s < S; s++) {
        const v = level - (marked && s % 4 === 2 ? 0.35 : 0) + (rng.float() - 0.5) * 0.1;
        value[i * S + s] = Math.round(Math.min(1, Math.max(0, v)) * 1000) / 1000;
      }
    });
    items.set(chr, {
      form: 'items',
      start: Int32Array.from(list, (g) => g.start),
      end: Int32Array.from(list, (g) => g.end),
      states: S,
      value,
    });
  }
  return levels(a, items, GeneMethylationBin);
}
