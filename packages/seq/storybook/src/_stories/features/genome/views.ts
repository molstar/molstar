/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { BedRuntime } from '@molstar/seqx-core/runtimes/bed';
import { registry as genomeBedRegistry } from '@molstar/seqx-examples/genomes/bed';
import { CytobandsFile } from '@molstar/seqx-examples/genomes/files';
import type { ExampleView } from '@molstar/seqx-examples';
import { blockSpan } from '@molstar/seqx-core/view';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { SeqViewer } from '@molstar/seqx-viewer/viewer';

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
  const p = b.data.assembly('peach', assembly('peach', peach, { order: 'source' }));
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

/**
 * One chromosome at a time (the story's viewer config): a pan, a zoom and a fit kept in the chromosome the window is in
 * (IGV's, UCSC's); the navigator's zoom strip under a window narrower than 200 px (a chromosome's at the stop, so the
 * strip shows the chromosome in view).
 */
export const OneAtATime = [
  [viewerItem('navigate.acrossBlocks'), 'stop'],
  [viewerItem('navigator.zoomStrip'), { minPx: 200, heightPx: 16 }],
] as const;

/**
 * A chromosome at a time: the genome's section with a navigator over its blocks (every chromosome and its label, a click
 * on one a jump there; under it the zoom strip of the chromosome holding the window); the viewer's
 * `navigate.acrossBlocks: stop` keeps a pan and a zoom in the chromosome. (An ideogram's bands over every chromosome:
 * the cytobands story.)
 */
function chromosomeAtATime(): Node {
  const b = createSeqSpec({ title: 'A chromosome at a time' });
  const g = b.data.assembly('toy', assembly(AssemblyName, toyKaryotype()));
  b.views.view('main', (v) => {
    const s = genome(v, 'genome', g);
    v.navigator('genome-nav', { of: s, axes: 'positions', content: 'axis', size: 28 } as never);
  });
  return b.build();
}

/** Where the square opens: its view on chr3 (a window off the whole, so a click on another cell centres there). */
export function openSquare(viewer: SeqViewer): void {
  const s = viewer.view.resolved.value?.sections[0];
  const b = s?.blocks?.[2];
  if (s && b) viewer.view.setRange(s.address, b.offset ?? 0, (b.offset ?? 0) + blockSpan(b));
}

/**
 * The square over blocks × blocks: one section on both axes (`of: [s, s]`), its navigator drawing the assembly's
 * chromosome pairs as a lattice, the cells the section shows plain and the rest washed, the window one rectangle.
 */
function squareOverBlocks(): Node {
  const b = createSeqSpec({ title: 'The square over blocks × blocks' });
  const g = b.data.assembly('toy', assembly(AssemblyName, toyKaryotype().slice(0, 12)));
  b.views.view('main', (v) => {
    const s = genome(v, 'genome', g);
    v.navigator('square', { of: [s, s], axes: 'both', system: g, size: [220, 220] } as never);
  });
  return b.build();
}

/**
 * A ring with `navigator.switch: include`: the genome's ring and a navigator over the assembly; a click on a chromosome
 * shows it alone on the ring (`section.include`, one undo step), a click on the navigator's title or Home all again.
 */
function ringSwitch(): Node {
  const b = createSeqSpec({ title: 'A ring that switches' });
  const g = b.data.assembly('toy', assembly(AssemblyName, toyKaryotype()));
  b.views.view('main', (v) => {
    const s = genome(v, 'ring', g, { layout: 'rings' });
    v.navigator('chromosomes', {
      of: s,
      axes: 'positions',
      content: 'axis',
      system: g,
      switch: 'include',
      size: 28,
    } as never);
  });
  return b.build();
}

/** The generated BED over the toy assembly: features on each chromosome, its names without the `chr` prefix. */
export function toyBed(seed = 3): string {
  const rng = new PCG(seed);
  const lines: string[] = [];
  for (const c of toyKaryotype()) {
    if (c.cyclic) continue;
    // a dozen features a chromosome, 0.5 to 3 Mb, on both strands
    for (let i = 0; i < 12; i++) {
      const start = Math.floor(((i + (rng.int() % 1000) / 1000) * c.length) / 12);
      const end = Math.min(c.length, start + 500_000 + (rng.int() % 2_500_000));
      lines.push([c.name.slice(3), start, end, `${c.name.slice(3)}.f${i + 1}`, 0, i % 3 ? '+' : '-'].join('\t'));
    }
  }
  return `${lines.join('\n')}\n`;
}

/**
 * A BED over the assembly's rows: one file over every chromosome (its names `1`, `X`, the assembly's `chr1`, `chrX`:
 * the names rule maps them), read by region a chromosome at a time, each block drawing its own features.
 */
function bedOverRows(): Node {
  const b = createSeqSpec({ title: 'A BED over the assembly' });
  const g = b.data.assembly('toy', assembly(AssemblyName, toyKaryotype()));
  const bed = b.data.source('features-bed', { format: 'bed', value: toyBed() as never });
  const features = b.data.ranges('features', { target: g, source: bed });
  b.views.view('main', (v) =>
    genome(v, 'genome', g, {}, (s) =>
      s.track('features', { label: 'Features' }).layer('features', { type: 'ranges', data: features }),
    ),
  );
  return b.build();
}

/** GRCh38's primary chromosomes and their lengths (the cytobands' last band ends there). */
export const Grch38: readonly { name: string; length: number; cyclic?: boolean }[] = [
  ...[
    248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
    50818468,
  ].map((length, i) => ({ name: `chr${i + 1}`, length })),
  { name: 'chrX', length: 156_040_895 },
  { name: 'chrY', length: 57_227_415 },
  { name: 'chrM', length: 16_569, cyclic: true },
];

/**
 * The cytobands over every chromosome: hg38's cytoBandIdeo (the recorded file) over GRCh38's primary assembly, read
 * by sequence, as a bands track on every block and as the ideogram navigator over the genome.
 */
function cytobandsOverGenome(): Node {
  const b = createSeqSpec({ title: 'Cytobands over the genome' });
  const g = b.data.assembly('hg38', assembly('GRCh38', Grch38, { aliases: ['hg38'] }));
  const src = b.data.source('cytobands', {
    format: 'bed',
    uri: CytobandsFile,
    fields: { stain: '4' },
    attribution: 'UCSC Genome Browser, hg38 cytoBandIdeo',
    licence: 'LicenseRef-UCSC-freely-usable',
  } as never);
  const bands = b.data.ranges('bands', { target: g, source: src });
  b.views.view('main', (v) =>
    genome(v, 'genome', g, { bands }, (s) =>
      s.track('bands', { label: 'Bands' }).layer('bands', { type: 'bands', data: bands }),
    ),
  );
  return b.build();
}

/**
 * Links written in the spec over the toy assembly, a chord ring: each row names its two sequences (`a_name`, `b_name`,
 * some without `chr`: the names rule maps them), a link between two chromosomes a chord between their arcs.
 */
export function inlineLinks(seed = 5) {
  const rng = new PCG(seed);
  const chromosomes = toyKaryotype().filter((c) => !c.cyclic);
  const rows = {
    a_start: [] as number[],
    a_end: [] as number[],
    b_start: [] as number[],
    b_end: [] as number[],
    a_name: [] as string[],
    b_name: [] as string[],
  };
  const name = (c: { name: string }, i: number) => (i % 3 ? c.name : c.name.slice(3));
  for (let i = 0; i < 60; i++) {
    const a = chromosomes[rng.int() % chromosomes.length];
    const b = chromosomes[(rng.int() % 5 ? rng.int() : 0) % chromosomes.length];
    const as = Math.floor((a.length * (rng.int() % 1000)) / 1000);
    const bs = Math.floor((b.length * (rng.int() % 1000)) / 1000);
    rows.a_start.push(Math.min(as, a.length - 1_000_000));
    rows.a_end.push(Math.min(as, a.length - 1_000_000) + 1_000_000);
    rows.b_start.push(Math.min(bs, b.length - 1_000_000));
    rows.b_end.push(Math.min(bs, b.length - 1_000_000) + 1_000_000);
    rows.a_name.push(name(a, i));
    rows.b_name.push(name(b, i + 1));
  }
  return rows;
}

function inlineLinksRing(): Node {
  const b = createSeqSpec({ title: 'Links written in the spec, on a ring' });
  const g = b.data.assembly('toy', assembly(AssemblyName, toyKaryotype()));
  const links = b.data.pairs('links', { a: g, rows: inlineLinks() });
  b.views.view('main', (v) =>
    genome(v, 'ring', g, { layout: 'rings' }, (s) =>
      s.track('links', { label: 'Links' }).layer('links', { type: 'pairs', data: links, legend: 'none' }),
    ),
  );
  return b.build();
}

/** Two short chromosomes (a toy: a position an item, so the items show at the fit). */
export const TwoChromosomes = [
  { name: 'chrI', length: 48 },
  { name: 'chrII', length: 32 },
];

/**
 * A matrix written in the spec over two chromosomes: four samples, a value an item and sample; `name` says each item's
 * chromosome (`I` is `chrI`), a chromosome's items its positions in order.
 */
function inlineMatrix(): Node {
  const b = createSeqSpec({ title: 'A matrix written in the spec, over two chromosomes' });
  const g = b.data.assembly('toy', assembly('Toy2', TwoChromosomes, { order: 'source' }));
  const rng = new PCG(11);
  const states = ['s1', 's2', 's3', 's4'];
  const values: number[] = [];
  const name: string[] = [];
  for (const c of TwoChromosomes)
    for (let i = 0; i < c.length; i++) {
      name.push(i % 2 ? c.name : c.name.slice(3));
      // a methylation-like level, smooth along the chromosome, a sample's offset
      for (let k = 0; k < states.length; k++)
        values.push(Math.max(0, Math.min(1, 0.5 + 0.4 * Math.sin(i / 6 + k) + ((rng.int() % 100) - 50) / 500)));
    }
  const m = b.data.matrix('levels', { target: g, values, states, name });
  b.views.view('main', (v) =>
    genome(v, 'genome', g, {}, (s) =>
      s.track('levels', { label: 'Levels', height: 80 }).layer('levels', { type: 'matrix', data: m }),
    ),
  );
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
  {
    dataset: 'genome',
    id: 'chromosome-at-a-time',
    title: 'A chromosome at a time',
    summary:
      'A navigator over the genome’s blocks (every chromosome and its name; the bands are the cytobands story’s): a click on one a jump there, the zoom strip of the chromosome in view under it; the readout names the chromosome and the position.',
    build: chromosomeAtATime,
  },
  {
    dataset: 'genome',
    id: 'square',
    title: 'The square over blocks × blocks',
    summary:
      'A navigator on both axes over the genome’s twelve chromosomes: their pairs as a lattice sized by length, named on both axes, the window one rectangle; a click centres there, a double-click fits a cell’s pair. Its cells’ reads are region reads’ (not drawn yet).',
    build: squareOverBlocks,
  },
  {
    dataset: 'genome',
    id: 'ring-switch',
    title: 'A ring that switches',
    summary:
      'The genome’s ring and a navigator over the assembly with `switch: include`: a click on a chromosome shows it alone on the ring (one undo step); its title or Home shows every one again.',
    build: ringSwitch,
  },
  {
    dataset: 'genome',
    id: 'bed-over-rows',
    title: 'A BED over the assembly',
    summary:
      'One BED file over every chromosome of the toy assembly (its names without `chr`: the names rule maps them), read by region a chromosome at a time: each block draws its own features.',
    build: bedOverRows,
    registry: () => [{ kind: 'format-runtime', runtime: BedRuntime }],
  },
  {
    dataset: 'genome',
    id: 'cytobands',
    title: 'Cytobands over the genome',
    summary:
      'hg38’s cytobands (UCSC’s cytoBandIdeo, the recorded file) over GRCh38’s 25 primary sequences, read by sequence: a bands track on every chromosome and the ideogram navigator over the genome; the file’s 686 other contigs are said once.',
    build: cytobandsOverGenome,
    registry: (o) => genomeBedRegistry(o),
  },
  {
    dataset: 'genome',
    id: 'inline-links-ring',
    title: 'Links written in the spec, on a ring',
    summary:
      'Sixty links written in the spec over the toy assembly, each row naming its two chromosomes (`a_name`, `b_name`; some without `chr`, the names rule maps them): a chord ring, a link between two chromosomes a chord between their arcs.',
    build: inlineLinksRing,
  },
  {
    dataset: 'genome',
    id: 'inline-matrix',
    title: 'A matrix written in the spec, over two chromosomes',
    summary:
      'Four samples’ levels written in the spec over two short chromosomes, `name` an item’s chromosome: each block draws its own items, a chromosome’s items its positions in order.',
    build: inlineMatrix,
  },
];
