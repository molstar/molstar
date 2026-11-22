/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import { categoryColumn } from '@molstar/core/repr/tree';
import type { Coverage, PairsList, SeqData, SeqPair, SystemInfo } from '@molstar/seqx-core/data';

// Three generated toy genomes and their alignments (Features/Across sections): peach, grape made from peach by an
// inversion, a translocation and a gap, and plum made from grape by a second translocation and a fusion. An alignment
// is its blocks (gapless runs of a few kb, a block a pair: its extent on each genome, the strand of its `b` end), served
// as pairs in the `seq` format over both genomes' chromosomes (their names distinct, as their assemblies name them:
// peach's `Pp01`, grape's `chr1`, plum's `Pd1`), each found from either end. Small (megabases, hundreds of blocks) and
// deterministic on core's `PCG`. Why not sided: a bridge reads a source's pairs a chromosome at a time without a side
// yet (the lane's next step), so the names must not collide.

export interface Chromosome {
  readonly name: string;
  readonly length: number;
}

/** A stretch of one genome carried to another: `[a0, a1)` on `a` to `b` from `b0`, read backwards when `strand` −1. */
export interface Segment {
  readonly a: string;
  readonly a0: number;
  readonly a1: number;
  readonly b: string;
  readonly b0: number;
  readonly strand: 1 | -1;
}

/** A genome: its assembly's name, its chromosomes. */
export interface ToyGenome {
  readonly name: string;
  readonly chromosomes: readonly Chromosome[];
}

export const Peach: ToyGenome = {
  name: 'Peach',
  chromosomes: [
    { name: 'Pp01', length: 3_000_000 },
    { name: 'Pp02', length: 2_600_000 },
    { name: 'Pp03', length: 2_200_000 },
    { name: 'Pp04', length: 1_800_000 },
  ],
};

/** Peach's `Pp02` read backwards in grape's `chr2` (an inversion). */
export const Inversion: readonly [number, number] = [800_000, 1_600_000];
/** Peach's `Pp01` stretch grape lost (a gap: no ribbon). */
export const Gap: readonly [number, number] = [1_000_000, 1_300_000];
/** Where peach's `Pp03` breaks: its tail is grape's `chr4` (a translocation). */
export const Translocation = 1_400_000;

export const Grape: ToyGenome = {
  name: 'Grape',
  chromosomes: [
    { name: 'chr1', length: 2_700_000 },
    { name: 'chr2', length: 2_600_000 },
    { name: 'chr3', length: 1_400_000 },
    { name: 'chr4', length: 1_400_000 },
    { name: 'chr5', length: 1_800_000 },
  ],
};

/** Peach to grape, stretch by stretch. */
export const PeachToGrape: readonly Segment[] = [
  { a: 'Pp01', a0: 0, a1: Gap[0], b: 'chr1', b0: 0, strand: 1 },
  { a: 'Pp01', a0: Gap[1], a1: 3_000_000, b: 'chr1', b0: Gap[0], strand: 1 },
  { a: 'Pp02', a0: 0, a1: Inversion[0], b: 'chr2', b0: 0, strand: 1 },
  { a: 'Pp02', a0: Inversion[0], a1: Inversion[1], b: 'chr2', b0: Inversion[0], strand: -1 },
  { a: 'Pp02', a0: Inversion[1], a1: 2_600_000, b: 'chr2', b0: Inversion[1], strand: 1 },
  { a: 'Pp03', a0: 0, a1: Translocation, b: 'chr3', b0: 0, strand: 1 },
  // Pp03's tail went to grape's chr4, before 600 kb of grape's own
  { a: 'Pp03', a0: Translocation, a1: 2_200_000, b: 'chr4', b0: 0, strand: 1 },
  { a: 'Pp04', a0: 0, a1: 1_800_000, b: 'chr5', b0: 0, strand: 1 },
];

export const Plum: ToyGenome = {
  name: 'Plum',
  chromosomes: [
    { name: 'Pd1', length: 1_900_000 },
    { name: 'Pd2', length: 2_600_000 },
    { name: 'Pd3', length: 2_200_000 },
    { name: 'Pd4', length: 3_200_000 },
  ],
};

/** Grape to plum: grape's `chr1` tail joins its `chr5` (plum's `Pd4`), its `chr3` and `chr4` fuse (plum's `Pd3`). */
export const GrapeToPlum: readonly Segment[] = [
  { a: 'chr1', a0: 0, a1: 1_900_000, b: 'Pd1', b0: 0, strand: 1 },
  { a: 'chr1', a0: 1_900_000, a1: 2_700_000, b: 'Pd4', b0: 2_400_000, strand: -1 },
  { a: 'chr2', a0: 0, a1: 2_600_000, b: 'Pd2', b0: 0, strand: 1 },
  { a: 'chr3', a0: 0, a1: 1_400_000, b: 'Pd3', b0: 0, strand: 1 },
  { a: 'chr4', a0: 0, a1: 800_000, b: 'Pd3', b0: 1_400_000, strand: 1 },
  { a: 'chr5', a0: 0, a1: 1_800_000, b: 'Pd4', b0: 0, strand: 1 },
];

/** Peach against itself: three duplicated stretches (one on its own chromosome, one elsewhere, one reversed). */
export const PeachDuplications: readonly Segment[] = [
  { a: 'Pp01', a0: 200_000, a1: 600_000, b: 'Pp03', b0: 1_700_000, strand: 1 },
  { a: 'Pp02', a0: 100_000, a1: 400_000, b: 'Pp04', b0: 1_200_000, strand: -1 },
  { a: 'Pp02', a0: 1_900_000, a1: 2_100_000, b: 'Pp02', b0: 2_300_000, strand: 1 },
];

/** An alignment's blocks: a pair a gapless run, `pitch` apart on average, of both genomes' extents. */
export interface Blocks {
  readonly a: string[];
  readonly aStart: number[];
  readonly aEnd: number[];
  readonly b: string[];
  readonly bStart: number[];
  readonly bEnd: number[];
  readonly strand: (1 | -1)[];
}

/**
 * The blocks of `segments`: each segment walked from its start by `pitch`, a block 35 to 75 % of it (the rest the
 * alignment's gap), its image on `b` by the segment (backwards on strand −1). Seeded.
 */
export function blocksOf(segments: readonly Segment[], pitch: number, seed: number): Blocks {
  const rng = new PCG(seed);
  const out: Blocks = { a: [], aStart: [], aEnd: [], b: [], bStart: [], bEnd: [], strand: [] };
  for (const s of segments)
    for (let at = s.a0; at + pitch <= s.a1; at += pitch) {
      const a0 = at + Math.floor(rng.float() * pitch * 0.15);
      const a1 = Math.min(s.a1, a0 + Math.floor(pitch * (0.35 + 0.4 * rng.float())));
      const b0 = s.strand > 0 ? s.b0 + (a0 - s.a0) : s.b0 + (s.a1 - a1);
      out.a.push(s.a);
      out.aStart.push(a0);
      out.aEnd.push(a1);
      out.b.push(s.b);
      out.bStart.push(b0);
      out.bEnd.push(b0 + (a1 - a0));
      out.strand.push(s.strand);
    }
  return out;
}

/** Blocks as a pairs table (`rows` of them; all when absent). */
function listOf(bl: Blocks, rows?: readonly number[]): PairsList {
  const at = rows ?? bl.a.map((_, i) => i);
  return {
    form: 'list',
    aSystem: categoryColumn(at.map((i) => bl.a[i])),
    bSystem: categoryColumn(at.map((i) => bl.b[i])),
    aStart: Int32Array.from(at, (i) => bl.aStart[i]),
    aEnd: Int32Array.from(at, (i) => bl.aEnd[i]),
    bStart: Int32Array.from(at, (i) => bl.bStart[i]),
    bEnd: Int32Array.from(at, (i) => bl.bEnd[i]),
    aStrand: new Int8Array(at.length).fill(1),
    bStrand: Int8Array.from(at, (i) => bl.strand[i]),
    row: Int32Array.from(at),
  };
}

/**
 * An alignment's blocks as a `seq` source, level 1 alone (the store's summaries give the pairs their levels): one list
 * of both genomes' chromosomes (one genome against itself: its own), each request answering the blocks with an end on
 * its chromosome, either way.
 */
export function alignmentSource(upper: ToyGenome, lower: ToyGenome | undefined, bl: Blocks): SeqData {
  const systems: SystemInfo[] = [...upper.chromosomes, ...(lower?.chromosomes ?? [])].map((c) => ({
    name: c.name,
    length: c.length,
  }));
  const coverage: Coverage = {
    systems,
    shapes: systems.map((s) => ({ system: s.name, shape: 'pairs' as const })),
    levels: [1],
  };
  const pairs: SeqPair[] = systems.map((s) => {
    const rows: number[] = [];
    for (let i = 0; i < bl.a.length; i++) if (bl.a[i] === s.name || bl.b[i] === s.name) rows.push(i);
    return {
      request: { system: s.name, start: 0, end: s.length!, level: 1, shapes: ['pairs'] },
      answer: { ok: true, tables: { pairs: listOf(bl, rows) } },
    };
  });
  return { coverage, pairs };
}

/** The blocks' pitch of the small alignments (a few hundred blocks: items at the fit) and of the dense one (bins). */
export const ItemsPitch = 24_000;
export const BinsPitch = 1_600;

/** Peach to grape, a few hundred blocks (with `pitch` 1,600: about 6,000). */
export const peachGrape = (pitch = ItemsPitch) => alignmentSource(Peach, Grape, blocksOf(PeachToGrape, pitch, 7));
/** Grape to plum. */
export const grapePlum = () => alignmentSource(Grape, Plum, blocksOf(GrapeToPlum, ItemsPitch, 11));
/** Peach against itself (its duplications). */
export const peachSelf = () => alignmentSource(Peach, undefined, blocksOf(PeachDuplications, 8_000, 5));
