/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import type { ExampleView } from '@molstar/seqx-examples';

// Toy genomes generated here, seeded (no downloads: the real dataset is the examples'), for the Genome coordinates
// stories: one assembly in a row and on a ring, two assemblies on one ring, one genome's chromosome in another's order
// with an inverted block reversed, and 3,000 contigs in runs.

/** A human-like toy karyotype: 22 autosomes from 250 Mb down to 46 Mb, X, Y and the mitochondrion, seeded. */
export function toyKaryotype(seed = 1): { name: string; length: number; cyclic?: boolean }[] {
  const rng = new PCG(seed);
  const out: { name: string; length: number; cyclic?: boolean }[] = [];
  for (let i = 1; i <= 22; i++)
    out.push({
      name: `chr${i}`,
      length: Math.round((250 - (i - 1) * 9.5) * 1e6 * (0.94 + 0.12 * ((rng.int() % 1000) / 1000))),
    });
  out.push({ name: 'chrX', length: 156_040_895 }, { name: 'chrY', length: 57_227_415 });
  out.push({ name: 'chrM', length: 16_569, cyclic: true });
  return out;
}

/** `n` unplaced contigs of 1 to 200 kb, seeded, named as a draft's. */
export function toyContigs(n: number, seed = 7): { name: string; length: number }[] {
  const rng = new PCG(seed);
  return Array.from({ length: n }, (_, i) => ({
    name: `chrUn_KI${270_300 + i}v1`,
    length: 1_000 + (rng.int() % 199_000),
  }));
}

const AssemblyName = 'Toy38';

function oneAssembly(layout: 'row' | 'ring'): Node {
  const b = createSeqSpec({ title: `One assembly ${layout === 'ring' ? 'on a ring' : 'in a row'}` });
  const g = b.data.assembly('toy', assembly(AssemblyName, toyKaryotype(), { aliases: ['toy38'] }));
  b.views.view('main', (v) => {
    if (layout === 'ring') {
      // the ring, and its linear detail following it: equal block lists follow by identity in units
      genome(v, 'ring', g, { layout: 'rings' });
      const detail = genome(v, 'genome', g);
      void detail;
    } else genome(v, 'genome', g);
  });
  const spec = b.build();
  if (layout === 'ring') follow(spec, 'genome', '/views/main/ring');
  return spec;
}

/** Sets a section's `follow` (the builder's `genome()` writes none). */
function follow(spec: Node, section: string, leader: string): void {
  const view = spec.children!.find((c) => c.kind === 'views')!.children![0];
  const s = view.children!.find((c) => c.id === section)! as { params?: Node['params'] };
  s.params = { ...s.params, follow: leader };
}

function twoOnARing(): Node {
  const b = createSeqSpec({ title: 'Two assemblies on one ring' });
  const rng = new PCG(3);
  const peach = Array.from({ length: 8 }, (_, i) => ({ name: `Pp0${i + 1}`, length: 18e6 + (rng.int() % 30_000_000) }));
  const grape = Array.from({ length: 19 }, (_, i) => ({
    name: `chr${i + 1}`,
    length: 18e6 + (rng.int() % 16_000_000),
  }));
  const p = b.data.assembly('peach', assembly('peach', peach, { order: 'source' } as never));
  const g = b.data.assembly('grape', assembly('grape', grape));
  b.views.view('main', (v) => genome(v, 'genomes', [p, g], { layout: 'rings' }));
  return b.build();
}

/** A's four chromosomes, and B's: B's chr2 carries an inversion of A's chr2 between 40 and 70 Mb. */
const ALengths = [120e6, 95e6, 80e6, 60e6];

function bInAsOrder(): Node {
  const b = createSeqSpec({ title: 'B in A’s order' });
  const a = b.data.assembly(
    'a',
    assembly(
      'A',
      ALengths.map((length, i) => ({ name: `chr${i + 1}`, length })),
    ),
  );
  const bb = b.data.assembly(
    'b',
    assembly(
      'B',
      ALengths.map((length, i) => ({ name: `chr${i + 1}`, length: length + i * 1e6 })),
    ),
  );
  const seq = (name: string) => `/data/b/(sequence :name ${name})`;
  b.views.view('main', (v) => {
    genome(v, 'a', a);
    // B's blocks in A's order, as the builder's colinearBlocks writes them: chr2 cut at the inversion's edges, its
    // middle block reversed; chr4 at four times its length (a scaled block, labelled ×4)
    v.section('b', {}, (s) => {
      s.block('b1', { of: seq('chr1'), label: 'chr1' } as never);
      s.block('b2a', { of: seq('chr2'), range: [0, 40e6], label: 'chr2' } as never);
      s.block('b2b', { of: seq('chr2'), range: [40e6, 70e6], label: 'chr2 inv', direction: 'reversed' } as never);
      s.block('b2c', { of: seq('chr2'), range: [70e6, 96e6], label: 'chr2' } as never);
      s.block('b3', { of: seq('chr3'), label: 'chr3' } as never);
      s.block('b4', { of: seq('chr4'), label: 'chr4', scale: 4 } as never);
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' } as never);
    });
    void bb;
  });
  return b.build();
}

function contigs(): Node {
  const b = createSeqSpec({ title: '3,000 contigs in runs' });
  const g = b.data.assembly('draft', assembly('Toy draft', [...toyKaryotype(2).slice(0, 24), ...toyContigs(2976)]));
  b.views.view('main', (v) => {
    genome(v, 'ring', g, { layout: 'rings' });
    genome(v, 'genome', g);
  });
  return b.build();
}

export const GenomeViews: readonly ExampleView[] = [
  {
    dataset: 'genome',
    id: 'row',
    title: 'One assembly in a row',
    summary:
      'A toy assembly of 25 sequences as blocks end to end, each its own ruler in Mb, the gaps sharing at most a quarter of the width (12 px each at most).',
    build: () => oneAssembly('row'),
  },
  {
    dataset: 'genome',
    id: 'ring',
    title: 'One assembly on a ring',
    summary:
      'The same assembly as arcs clockwise from the top, labels outside, the hole titled with its name and length; the linear detail under it follows the ring’s window (equal block lists).',
    build: () => oneAssembly('ring'),
  },
  {
    dataset: 'genome',
    id: 'two-on-a-ring',
    title: 'Two assemblies on one ring',
    summary:
      'peach and grape on one ring: a group each, a 4° gap between them, an arc outside and the assembly’s name.',
    build: twoOnARing,
  },
  {
    dataset: 'genome',
    id: 'b-in-a-order',
    title: 'B in A’s order',
    summary:
      'Genome A, and B’s chromosomes in A’s order: chr2 cut at an inversion’s edges, its middle block reversed (its ruler counts down, its label ←); chr4 at four times its length (×4).',
    build: bInAsOrder,
  },
  {
    dataset: 'genome',
    id: 'contigs',
    title: '3,000 contigs in runs',
    summary:
      'A draft of 24 chromosomes and 2,976 contigs on a ring and in a row: the contigs under 2 px a run, its count where it has room; a contig a device pixel wide is still read and drawn.',
    build: contigs,
  },
];
