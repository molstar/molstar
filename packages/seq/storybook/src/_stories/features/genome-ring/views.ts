/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec } from '@molstar/seqx-core/lang';
import type { ExampleView } from '@molstar/seqx-examples';
import { RingViews, registry, ringSpec } from '@molstar/seqx-examples/genomes/ring/views';
import { rondo, type RondoLay, rondoRing } from '@molstar/seqx-ext-rondo';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';
import { toyContigs, toyKaryotype } from '../genome/views.js';

// The genome ring's stories (Design › Genome ring, its defaults): the showcase's views (`RingViews`) and the
// defaults its data does not hold, each through `rondoRing` and drawn by a viewer composed with `rondoViewer`: two
// assemblies on one ring (synteny), a crowded draft (3,000 contigs), the links pending, failed and absent.

export const RingSectionId = 'ring';

const ring = (id: string) => RingViews.find((v) => v.id === id)!;

/** The ribbons' lay of a viewer's last frame (the first `rondo:ribbons` layer). */
export function ribbonsOf(viewer: SeqViewer): RondoLay | undefined {
  for (const s of viewer.frame?.layout.sections ?? [])
    for (const l of s.sectionLays) if (l.type.name === 'rondo:ribbons') return l.lay as RondoLay;
  return undefined;
}

/** The strongest placed stem's mouth, two px inward: where a pointer hovers it (default 6). */
export function strongestStem(l: RondoLay): { readonly stem: number; readonly x: number; readonly y: number } {
  let best = -1;
  for (let s = 0; s < l.ns; s++) if (l.stemPlaced[s] && (best < 0 || l.stemCount[s] > l.stemCount[best])) best = s;
  l.frame.point((l.stemT0[best] + l.stemT1[best]) / 2, 2);
  return { stem: best, x: l.frame.px, y: l.frame.py };
}

/** The story's opening for default 6: once the ribbons are laid out, the pointer on the strongest stem. */
export function hoverStrongestStem(viewer: SeqViewer): void {
  const l = ribbonsOf(viewer);
  if (!l || l.ns === 0) return;
  const { x, y } = strongestStem(l);
  viewer.input.move(x, y);
  viewer.draw();
}

/** Two assemblies on one ring (EC-210, peach and grape): synteny blocks between them as held pairs, one inverted run. */
function twoOnARing(): Node {
  const b = createSeqSpec({ title: 'Two assemblies on one ring', extensions: [rondo] });
  const rng = new PCG(3);
  const peach = Array.from({ length: 8 }, (_, i) => ({ name: `Pp0${i + 1}`, length: 18e6 + (rng.int() % 30_000_000) }));
  const grape = Array.from({ length: 19 }, (_, i) => ({
    name: `chr${i + 1}`,
    length: 18e6 + (rng.int() % 16_000_000),
  }));
  const p = b.data.assembly('peach', assembly('peach', peach, { order: 'source' }));
  const g = b.data.assembly('grape', assembly('grape', grape));
  const rows = {
    a_name: [] as string[],
    b_name: [] as string[],
    a_start: [] as number[],
    a_end: [] as number[],
    b_start: [] as number[],
    b_end: [] as number[],
  };
  // each peach chromosome's blocks on two or three grape chromosomes in order; Pp03's middle inverted on chr7
  for (let i = 0; i < peach.length; i++) {
    const partners = [(i * 2) % grape.length, (i * 2 + 1) % grape.length, (i * 5 + 3) % grape.length];
    const n = 40;
    for (let k = 0; k < n; k++) {
      const j = partners[k < 18 ? 0 : k < 32 ? 1 : 2];
      const at = Math.floor(((k + 0.5) / n) * (peach[i].length - 400_000));
      const share = (k % 18) / 18;
      const inverted = i === 2 && k >= 18 && k < 32;
      const bt = Math.floor((inverted ? 1 - share : share) * (grape[j].length - 400_000));
      rows.a_name.push(peach[i].name);
      rows.b_name.push(grape[j].name);
      rows.a_start.push(at);
      rows.a_end.push(at + 300_000 + (rng.int() % 100_000));
      rows.b_start.push(bt);
      rows.b_end.push(bt + 300_000 + (rng.int() % 100_000));
    }
  }
  const links = b.data.pairs('synteny', { a: p, b: g, rows } as never);
  b.views.view('main', (v) =>
    v
      .section(RingSectionId, { of: [p, g], layout: 'rings', center: 'none' } as never)
      .track('links', { label: 'Synteny', height: 4 })
      .layer('ribbons', { type: 'rondo:ribbons', fallback: 'pairs', data: links } as never),
  );
  return b.build();
}

/** A draft of 24 chromosomes and 2,976 contigs on the ring: links between chromosomes, and a few to contigs. */
function crowded(): Node {
  const b = createSeqSpec({ title: '3,000 contigs on the ring', extensions: [rondo] });
  const chromosomes = toyKaryotype(2).slice(0, 24);
  const contigs = toyContigs(2976);
  const g = b.data.assembly('draft', assembly('Toy draft', [...chromosomes, ...contigs]));
  const rng = new PCG(11);
  const rows = {
    a_name: [] as string[],
    b_name: [] as string[],
    a_start: [] as number[],
    a_end: [] as number[],
    b_start: [] as number[],
    b_end: [] as number[],
  };
  const end = (c: { length: number }) => Math.min(c.length - 1, 5000);
  for (let i = 0; i < 3000; i++) {
    const a = chromosomes[rng.int() % chromosomes.length];
    // one in five to a contig (a scaffold's join), the rest between chromosomes
    const c = rng.int() % 5 === 0 ? contigs[rng.int() % contigs.length] : chromosomes[rng.int() % chromosomes.length];
    const as = rng.int() % (a.length - 5000);
    const bs = rng.int() % Math.max(1, c.length - 5000);
    rows.a_name.push(a.name);
    rows.b_name.push(c.name);
    rows.a_start.push(as);
    rows.a_end.push(as + 5000);
    rows.b_start.push(bs);
    rows.b_end.push(bs + end(c));
  }
  const links = b.data.pairs('links', { a: g, rows } as never);
  b.views.view('main', (v) => {
    rondoRing(v, { id: RingSectionId, assembly: g, links, square: false });
  });
  return b.build();
}

/**
 * The stories' views, a default each (the design's numbers): the showcase's ring, chromosome and pair, the same ring
 * without the extension, the defaults the showcase's data does not hold.
 */
export const RingStoryViews: readonly ExampleView[] = [
  ring('ring'),
  ring('ring-chromosome'),
  ring('ring-pair'),
  {
    ...ring('ring'),
    id: 'ring-two-assemblies',
    title: 'Two assemblies on one ring',
    summary:
      'Peach’s 8 chromosomes and grape’s 19 on one ring, a group gap between them; synteny blocks as ribbons from a peach stem to a grape stem, Pp03’s middle inverted on its partner (the ribbon crossed).',
    build: twoOnARing,
  },
  {
    ...ring('ring'),
    id: 'ring-crowded',
    title: '3,000 contigs on the ring',
    summary:
      'A draft’s 24 chromosomes and 2,976 contigs: the contigs under 2 px a run; stems only on blocks a device px wide; the links to contigs end at their places in the run.',
    build: crowded,
  },
  {
    ...ring('ring'),
    id: 'ring-pending',
    title: 'The links pending',
    summary:
      'The links’ answer never comes (an hour’s delay): the bands drawn, the hole empty with the pending hatch, no stems.',
    registry: (o = {}) => registry({ ...o, delayMs: 3_600_000 }),
  },
  {
    ...ring('ring'),
    id: 'ring-failed',
    title: 'The links failed',
    summary:
      'The links’ file is not there: its issue said once; no ribbon, no hatch waiting (a source that never opened is failed), its ARIA “links failed”.',
    build: () => ringSpec({ failed: true, title: 'The links failed' }),
  },
  {
    ...ring('ring'),
    id: 'ring-no-links',
    title: 'No links',
    summary: 'No pairs at all: an empty hole, its ARIA description “no links”.',
    build: () => ringSpec({ links: 0, title: 'No links' }),
  },
];

export const ringStoryView = (id: string): ExampleView => RingStoryViews.find((v) => v.id === id)!;
