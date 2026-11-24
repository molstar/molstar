/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn, type Node } from '@molstar/core/repr/tree';
import { binPairs, type PairsList, type SeqData, type SeqPair } from '@molstar/seqx-core/data';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import type { Chromosome, Segment, ToyGenome } from './genomes.js';

// A liftover (Features/Across sections; EC-212): a toy assembly and its next version (`Rose v1`, `Rose v2`, their
// chromosomes named alike), the chain between them generated as gapless runs (a duplication, a deletion, an inversion
// and a translocation), served as sided pairs in the `seq` format (side `a` the old assembly, `b` the new), and the old
// assembly's genes. The figure: the old assembly above with its genes, the chain's runs as ribbons, the new assembly
// below with the old genes drawn through the chain. Deterministic on core's `PCG`; the truths (each gene's images) are
// computed here, from the runs, so the tests read them.

export const RoseV1: ToyGenome = {
  name: 'Rose v1',
  chromosomes: [
    { name: 'chr1', length: 3_000_000 },
    { name: 'chr2', length: 2_400_000 },
    { name: 'chr3', length: 1_800_000 },
  ],
};

export const RoseV2: ToyGenome = {
  name: 'Rose v2',
  chromosomes: [
    { name: 'chr1', length: 3_600_000 },
    { name: 'chr2', length: 2_100_000 },
    { name: 'chr3', length: 2_400_000 },
  ],
};

/** v1's `chr1` stretch v2 holds twice: in place, and again on `chr3` (a duplication). */
export const Duplicated: readonly [number, number] = [1_200_000, 2_000_000];
/** v1's `chr2` stretch v2 lost (a deletion: no run). */
export const Deleted: readonly [number, number] = [600_000, 900_000];
/** v1's `chr2` stretch v2 holds reversed (an inversion). */
export const Inverted: readonly [number, number] = [900_000, 1_700_000];
/** v1's `chr3` tail, at the end of v2's `chr1` (a translocation). */
export const Translocated: readonly [number, number] = [1_200_000, 1_800_000];

/** v1 to v2, stretch by stretch (`a` v1's names, `b` v2's). */
export const RoseChain: readonly Segment[] = [
  { a: 'chr1', a0: 0, a1: 3_000_000, b: 'chr1', b0: 0, strand: 1 },
  { a: 'chr1', a0: Duplicated[0], a1: Duplicated[1], b: 'chr3', b0: 1_600_000, strand: 1 },
  { a: 'chr2', a0: 0, a1: Deleted[0], b: 'chr2', b0: 0, strand: 1 },
  { a: 'chr2', a0: Inverted[0], a1: Inverted[1], b: 'chr2', b0: Deleted[0], strand: -1 },
  { a: 'chr2', a0: Inverted[1], a1: 2_400_000, b: 'chr2', b0: 1_400_000, strand: 1 },
  { a: 'chr3', a0: 0, a1: Translocated[0], b: 'chr3', b0: 0, strand: 1 },
  { a: 'chr3', a0: Translocated[0], a1: Translocated[1], b: 'chr1', b0: 3_000_000, strand: 1 },
];

/** A run of the chain: `[from, from + length)` on v1's `a` to v2's `b` from `to`, read backwards on strand −1. */
export interface ChainRun {
  readonly a: string;
  readonly from: number;
  readonly b: string;
  readonly to: number;
  readonly length: number;
  readonly strand: 1 | -1;
}

/**
 * The chain's runs: each segment cut into gapless runs of 60 to 140 kb with gaps of up to 2 kb between them (an
 * unaligned stretch: a gap in the ribbons), carried as the segment carries them. Seeded.
 */
export function chainRuns(segments: readonly Segment[] = RoseChain, seed = 3): ChainRun[] {
  const rng = new PCG(seed);
  const out: ChainRun[] = [];
  for (const s of segments) {
    let at = s.a0;
    while (at < s.a1) {
      const length = Math.min(s.a1 - at, 60_000 + Math.floor(rng.float() * 80_000));
      const to = s.strand > 0 ? s.b0 + (at - s.a0) : s.b0 + (s.a1 - at - length);
      out.push({ a: s.a, from: at, b: s.b, to, length, strand: s.strand });
      at += length + Math.floor(rng.float() * 2_000);
    }
  }
  return out;
}

/** A gene of v1: its id, chromosome, extent, strand and coding (thick) part. */
export interface Gene {
  readonly id: string;
  readonly chr: string;
  readonly start: number;
  readonly end: number;
  readonly strand: '+' | '-';
  readonly thickStart: number;
  readonly thickEnd: number;
}

/** A gene's strand and its thick part, the middle four fifths (drawn thick; lifted, carried through the runs). */
const coding = (id: string, chr: string, start: number, end: number, strand: '+' | '-'): Gene => {
  const pad = Math.floor((end - start) / 10);
  return { id, chr, start, end, strand, thickStart: start + pad, thickEnd: end - pad };
};

/**
 * v1's genes: one in each event (`dup` in the duplicated stretch, `lost` in the deleted one, `inv` in the inversion,
 * `moved` in the translocated tail), then 36 more of 10 to 50 kb along the chromosomes. Seeded.
 */
export function roseGenes(seed = 9): Gene[] {
  const rng = new PCG(seed);
  const genes: Gene[] = [
    coding('dup', 'chr1', 1_500_000, 1_540_000, '+'),
    coding('lost', 'chr2', 700_000, 740_000, '-'),
    coding('inv', 'chr2', 1_200_000, 1_250_000, '+'),
    coding('moved', 'chr3', 1_500_000, 1_530_000, '-'),
  ];
  for (let i = 0; i < 36; i++) {
    const c = RoseV1.chromosomes[i % RoseV1.chromosomes.length];
    const length = 10_000 + Math.floor(rng.float() * 40_000);
    const start = Math.floor(rng.float() * (c.length - length));
    genes.push(coding(`g${String(i + 1).padStart(2, '0')}`, c.name, start, start + length, i % 3 ? '+' : '-'));
  }
  return genes;
}

/** An image on v2: its chromosome and extent. */
export interface Image {
  readonly chr: string;
  readonly start: number;
  readonly end: number;
}

/** The images of v1's `[start, end)` on `chr` through `runs`: each run's part, sorted, the adjacent joined. */
export function liftImages(runs: readonly ChainRun[], chr: string, start: number, end: number): Image[] {
  const out: Image[] = [];
  for (const r of runs) {
    if (r.a !== chr) continue;
    const lo = Math.max(start, r.from);
    const hi = Math.min(end, r.from + r.length);
    if (hi <= lo) continue;
    if (r.strand > 0) out.push({ chr: r.b, start: r.to + lo - r.from, end: r.to + hi - r.from });
    else out.push({ chr: r.b, start: r.to + r.length - (hi - r.from), end: r.to + r.length - (lo - r.from) });
  }
  out.sort((x, y) => x.chr.localeCompare(y.chr) || x.start - y.start);
  const joined: Image[] = [];
  for (const m of out) {
    const last = joined[joined.length - 1];
    if (last && last.chr === m.chr && last.end >= m.start)
      joined[joined.length - 1] = { ...last, end: Math.max(last.end, m.end) };
    else joined.push(m);
  }
  return joined;
}

/** The runs as a pairs list: `a` on v1, `b` on v2, a row each. */
function runsList(runs: readonly ChainRun[]): PairsList {
  return {
    form: 'list',
    aSystem: categoryColumn(runs.map((r) => r.a)),
    bSystem: categoryColumn(runs.map((r) => r.b)),
    aStart: Int32Array.from(runs, (r) => r.from),
    aEnd: Int32Array.from(runs, (r) => r.from + r.length),
    bStart: Int32Array.from(runs, (r) => r.to),
    bEnd: Int32Array.from(runs, (r) => r.to + r.length),
    aStrand: new Int8Array(runs.length).fill(1),
    bStrand: Int8Array.from(runs, (r) => r.strand),
    row: Int32Array.from(runs, (_, i) => i),
  };
}

/** The chain's levels: its runs at 1, bin pairs above. */
const Levels = [1, 16_384];

/**
 * The chain as a sided `seq` source: each chromosome listed once a side (`chr1` on both), a request on a side answered
 * with the runs whose end on that side is on its chromosome (the runtime cuts them to the region, from either end).
 */
export function chainSource(runs: readonly ChainRun[] = chainRuns()): SeqData {
  const list = runsList(runs);
  const sides = [
    ['a', RoseV1.chromosomes],
    ['b', RoseV2.chromosomes],
  ] as const;
  const pairs: SeqPair[] = [];
  for (const [side, chrs] of sides)
    for (const c of chrs)
      for (const level of Levels)
        pairs.push({
          request: { system: c.name, start: 0, end: c.length, level, shapes: ['pairs'], side },
          answer: { ok: true, tables: { pairs: level === 1 ? list : binPairs(list, level) } },
        });
  return {
    coverage: {
      systems: sides.flatMap(([side, chrs]) => chrs.map((c: Chromosome) => ({ ...c, side }))),
      shapes: RoseV1.chromosomes.map((c) => ({ system: c.name, shape: 'pairs' as const })),
      levels: Levels,
    },
    pairs,
  };
}

/** The sections' addresses (what the tests read). */
export const OldSection = '/views/main/v1';
export const NewSection = '/views/main/v2';
export const ChainBridge = '/views/main/chain';

/**
 * The figure: v1 above with its genes, the chain's runs as a bridge's ribbons (`map.source`: the map's runs are the
 * source's rows), v2 below with v1's genes lifted through the same map. The chain a source of `format` (`seq-across`),
 * named `generated:liftover`.
 */
export function liftover(format: string): Node {
  const b = createSeqSpec({ title: 'Rose v1 lifted to v2' });
  const v1 = b.data.assembly('v1', assembly(RoseV1.name, RoseV1.chromosomes));
  const v2 = b.data.assembly('v2', assembly(RoseV2.name, RoseV2.chromosomes));
  const chain = b.data.source('chain', { format, uri: 'generated:liftover' } as never);
  const map = b.data.map('v1-v2', { from: v1, to: v2, method: 'segments', source: chain } as never);
  const genes = roseGenes();
  const data = b.data.ranges('genes', {
    target: v1,
    rows: {
      id: genes.map((g) => g.id),
      target: genes.map((g) => g.chr),
      start: genes.map((g) => g.start),
      end: genes.map((g) => g.end),
      strand: genes.map((g) => g.strand),
      thick_start: Int32Array.from(genes, (g) => g.thickStart),
      thick_end: Int32Array.from(genes, (g) => g.thickEnd),
    },
  } as never);
  b.views.view('main', (v) => {
    v.section('v1', { of: v1 } as never, (s) => {
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' } as never);
      s.track('genes', { label: 'Genes (v1)' }).layer('genes', {
        type: 'ranges',
        data,
        line: 'gene',
        legend: 'none',
      } as never);
    });
    v.bridge('chain', { a: OldSection, b: NewSection, label: 'v1 to v2' } as never, (br) =>
      br.layer('runs', { type: 'pairs', data: map, legend: 'none' } as never),
    );
    v.section('v2', { of: v2 } as never, (s) => {
      s.track('genes', { label: 'Genes (lifted from v1)' }).layer('genes', {
        type: 'ranges',
        data,
        line: 'gene',
        legend: 'none',
      } as never);
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' } as never);
    });
  });
  return b.build();
}
