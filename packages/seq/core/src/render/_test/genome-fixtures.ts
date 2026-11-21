/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Genome fixtures (seeded, generated: no downloads): blocks of mixed direction and scale over sequences with letters
// and features, a toy assembly of chromosomes and many small contigs (runs), two assemblies on one section.

import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { seeded } from './fixtures.js';

/** Seeded DNA letters. */
export function dna(n: number, seed: number): string {
  const r = seeded(seed);
  let s = '';
  for (let i = 0; i < n; i++) s += 'ACGT'[Math.floor(r() * 4)];
  return s;
}

/** A block as `mixedSpec` writes it: its sequence, range, direction and scale. */
export interface MixedBlock {
  readonly seq: number;
  readonly range: readonly [number, number];
  readonly direction?: 'section' | 'reversed';
  readonly scale?: number;
}

/** The default mixed blocks: forward, reversed, scaled, reversed and scaled, on two sequences of 400 and 300. */
export const MixedBlocks: readonly MixedBlock[] = [
  { seq: 0, range: [0, 150] },
  { seq: 0, range: [150, 400], direction: 'reversed' },
  { seq: 1, range: [0, 100], scale: 3 },
  { seq: 1, range: [100, 300], direction: 'reversed', scale: 0.5 },
];

/** Its features: a feature a block's interior on sequence 0 (forward and reversed) and on sequence 1. */
export const MixedFeatures = [
  { seq: 0, start: 20, end: 60, id: 'f0' },
  { seq: 0, start: 200, end: 260, id: 'f1' },
  { seq: 1, start: 30, end: 70, id: 'f2' },
  { seq: 1, start: 150, end: 250, id: 'f3' },
];

/** Blocks of mixed direction and scale over two sequences, a letters, a features and a ruler track. */
export function mixedSpec(blocks: readonly MixedBlock[] = MixedBlocks, layout?: 'rings'): Node {
  const b = createSeqSpec({ title: 'mixed blocks' });
  const seqs = [
    b.data.sequence('s0', { alphabet: 'dna', letters: dna(400, 3), name: 'chrA' }),
    b.data.sequence('s1', { alphabet: 'dna', letters: dna(300, 5), name: 'chrB' }),
  ];
  const feats = seqs.map((s, k) => {
    const rows = MixedFeatures.filter((f) => f.seq === k);
    return b.data.ranges(`f${k}`, {
      target: s,
      rows: { start: rows.map((f) => f.start), end: rows.map((f) => f.end), id: rows.map((f) => f.id) },
    });
  });
  b.views.view('v', (v) => {
    const s = v.section('g', layout ? { layout } : {});
    blocks.forEach((x, i) =>
      s.block(`b${i}`, {
        of: seqs[x.seq],
        range: [...x.range],
        ...(x.direction ? { direction: x.direction } : {}),
        ...(x.scale !== undefined ? { scale: x.scale } : {}),
      } as never),
    );
    s.track('letters').layer('letters', { type: 'letters' });
    const t = s.track('features');
    feats.forEach((f, k) => t.layer(`f${k}`, { type: 'ranges', data: f }));
    s.track('ruler').layer('ruler', { type: 'ruler' });
  });
  return b.build();
}

/** The sequences of a toy genome: `chromosomes` of 2 to 250 Mb, then `contigs` of 1 to 200 kb, seeded. */
export function toyGenome(chromosomes: number, contigs: number, seed = 1): { name: string; length: number }[] {
  const r = seeded(seed);
  const out: { name: string; length: number }[] = [];
  for (let i = 0; i < chromosomes; i++)
    out.push({
      name: `chr${i + 1}`,
      length: Math.round(2e6 + (248e6 - 2e6) * (1 - i / Math.max(1, chromosomes)) * (0.7 + 0.3 * r())),
    });
  for (let i = 0; i < contigs; i++)
    out.push({ name: `chrUn_KI${270000 + i}v1`, length: Math.round(1e3 + 199e3 * r()) });
  return out;
}

/** A toy genome's section: one assembly in a row (or a ring), its ruler; a reversed or scaled block by name. */
export function genomeSpec(
  sequences: readonly { name: string; length: number }[],
  options: {
    layout?: 'rings';
    reversed?: readonly string[];
    scaled?: Readonly<Record<string, number>>;
    name?: string;
    tracks?: number;
  } = {},
): Node {
  const b = createSeqSpec({ title: 'genome' });
  b.data.assembly('g', assembly(options.name ?? 'toy', sequences, { order: 'source' } as never));
  b.views.view('v', (v) =>
    genome(v, 'genome', '/data/g', options.layout ? { layout: options.layout } : {}, (s) => {
      for (let k = 1; k < (options.tracks ?? 1); k++) s.track(`ruler${k}`).layer(`ruler${k}`, { type: 'ruler' });
      for (const name of options.reversed ?? [])
        s.block(`r_${name}`, { of: `/data/g/(sequence :name ${name})`, direction: 'reversed' } as never);
      for (const [name, scale] of Object.entries(options.scaled ?? {}))
        s.block(`x_${name}`, { of: `/data/g/(sequence :name ${name})`, scale } as never);
    }),
  );
  return b.build();
}

/** Two toy assemblies on one section (EC-210): `peach` of 8 and `grape` of 19 chromosomes. */
export function twoGenomesSpec(layout?: 'rings'): Node {
  const b = createSeqSpec({ title: 'two genomes' });
  const peach = toyGenome(8, 0, 2).map((s, i) => ({ name: `Pp0${i + 1}`, length: Math.round(s.length / 5) }));
  const grape = toyGenome(19, 0, 3).map((s) => ({ ...s, length: Math.round(s.length / 6) }));
  b.data.assembly('peach', assembly('peach', peach, { order: 'source' } as never));
  b.data.assembly('grape', assembly('grape', grape));
  b.views.view('v', (v) => genome(v, 'two', ['/data/peach', '/data/grape'], layout ? { layout } : {}));
  return b.build();
}
