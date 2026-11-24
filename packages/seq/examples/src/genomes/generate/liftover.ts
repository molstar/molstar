/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type ChainRun, liftImages } from '../../across/liftover.js';
import type { ToyAssembly } from './assemblies.js';

// Step 2's liftover (EC-211, EC-212): toy A′, a revision of toy A named alike, with a tandem
// duplication on chr1 (10 kb held twice), an insertion on chr2 (5 kb of new letters) and a deletion on chr3 (4 kb
// lost); the map between them as gapless runs (written as sided pairs in the `seq` format, side `a` A, `b` A′: their
// names collide; at level 1 and as bin pairs of 16,384, the level a bridge reads at the fit), and positions of A with
// their images on A′ computed from the runs (the truths the tests read).

/** A's chr1 stretch A′ holds twice, the copy right after it. */
export const Duplicated: readonly [number, number] = [200_000, 210_000];
/** Where A′'s chr2 has 5 kb A has not. */
export const Inserted: readonly [at: number, length: number] = [50_000, 5_000];
/** A's chr3 stretch A′ lost. */
export const Deleted: readonly [number, number] = [30_000, 34_000];

/** Toy A′: A's sequences, chr1 10 kb longer, chr2 5 kb longer, chr3 4 kb shorter. */
export function toyA2(a: ToyAssembly): ToyAssembly {
  const grow: Readonly<Record<string, number>> = {
    chr1: Duplicated[1] - Duplicated[0],
    chr2: Inserted[1],
    chr3: -(Deleted[1] - Deleted[0]),
  };
  return { name: 'Toy A′', sequences: a.sequences.map((s) => ({ ...s, length: s.length + (grow[s.name] ?? 0) })) };
}

/** The runs from A to A′ (a run a gapless stretch; `a` and `b` its sequence on each, named alike). */
export function liftoverRuns(a: ToyAssembly): ChainRun[] {
  const len = (n: string) => a.sequences.find((s) => s.name === n)!.length;
  const run = (c: string, from: number, to: number, length: number): ChainRun => ({
    a: c,
    from,
    b: c,
    to,
    length,
    strand: 1,
  });
  const [d0, d1] = Duplicated;
  const [i0, il] = Inserted;
  const [x0, x1] = Deleted;
  const out: ChainRun[] = [
    run('chr1', 0, 0, d1),
    run('chr1', d0, d1, d1 - d0),
    run('chr1', d1, d1 + (d1 - d0), len('chr1') - d1),
    run('chr2', 0, 0, i0),
    run('chr2', i0, i0 + il, len('chr2') - i0),
    run('chr3', 0, 0, x0),
    run('chr3', x1, x0, len('chr3') - x1),
  ];
  for (const s of a.sequences) if (!/^chr[123]$/.test(s.name)) out.push(run(s.name, 0, 0, s.length));
  return out;
}

/** The runs as the sided pairs' rows (`row` a run each). */
export function liftoverBlocks(runs: readonly ChainRun[]) {
  return runs.map((r, row) => ({
    a: r.a,
    aStart: r.from,
    aEnd: r.from + r.length,
    b: r.b,
    bStart: r.to,
    bEnd: r.to + r.length,
    strand: r.strand,
    row,
  }));
}

/** A probe of the liftover: a position of A and its images on A′ (none: lost; two: duplicated). */
export interface LiftoverTruth {
  readonly chr: string;
  readonly position: number;
  readonly images: readonly number[];
}

/** Probes in each event and beside them, each with its images through the runs. */
export function liftoverTruths(runs: readonly ChainRun[]): LiftoverTruth[] {
  const probes: [string, number][] = [
    ['chr1', 100_000],
    ['chr1', 205_000],
    ['chr1', 250_000],
    ['chr2', 40_000],
    ['chr2', 60_000],
    ['chr3', 20_000],
    ['chr3', 32_000],
    ['chr3', 50_000],
  ];
  return probes.map(([chr, position]) => ({
    chr,
    position,
    images: liftImages(runs, chr, position, position + 1).map((m) => m.start),
  }));
}
