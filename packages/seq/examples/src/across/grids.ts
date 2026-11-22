/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PCG } from '@molstar/core/data/util/hash-functions';
import type { Node } from '@molstar/core/repr/tree';
import { assembly, createSeqSpec, genome } from '@molstar/seqx-core/lang';
import { Grape, Peach } from './genomes.js';

// Grids across sections (Features/Across sections): a grid over two sections of blocks (a dot plot of peach against
// grape, contacts between two chromosomes' sections, EC-209) and a grid on categorical axes (a kinship grid over a
// genotype matrix's samples, EC-35). Generated, small, deterministic on core's `PCG`; the truths each story shows are
// computed here, so the tests read them.

/** The dot plot's section (x follows peach, y grape). */
export const DotsSection = '/views/main/dots';

/**
 * Peach and grape as genome sections, a grid following both drawing their alignment's blocks as dots and diagonals; the
 * alignment a source of `format` (`seq-across`), its blocks `generated:peach-grape`.
 */
export function dotPlot(format: string): Node {
  const b = createSeqSpec({ title: 'Peach against grape, a dot plot' });
  const peach = b.data.assembly('peach', assembly(Peach.name, Peach.chromosomes));
  const grape = b.data.assembly('grape', assembly(Grape.name, Grape.chromosomes));
  const src = b.data.source('alignment', { format, uri: 'generated:peach-grape' } as never);
  const aln = b.data.pairs('aln', { a: peach, b: grape, source: src });
  b.views.view('main', (v) => {
    genome(v, 'peach', peach);
    genome(v, 'grape', grape);
    v.section('dots', { layout: 'grid', follow: '/views/main/peach', follow_y: '/views/main/grape' } as never)
      .track('alignment', { label: 'Alignment' })
      .layer('aln', { type: 'pairs', data: aln, legend: 'none' } as never);
  });
  return b.build();
}

/** Two chromosomes' contact bins (a bin a position, 10 kb each): their lengths, and the hotspots between them. */
export const ContactBins = { a: 240, b: 180 } as const;
/** Where the two chromosomes touch: a bin pair each, the contacts' strongest. */
export const Hotspots: readonly (readonly [number, number])[] = [
  [60, 40],
  [170, 130],
];

/** The two chromosomes' inter-chromosomal contacts: a low noisy floor and a blob round each hotspot. Seeded. */
export function contactValues(): Float32Array {
  const r = new PCG(17);
  const { a, b } = ContactBins;
  const v = new Float32Array(a * b);
  for (let i = 0; i < a; i++)
    for (let j = 0; j < b; j++) {
      let x = 0.04 + 0.06 * r.float();
      for (const [hi, hj] of Hotspots) {
        const d2 = (i - hi) ** 2 + (j - hj) ** 2;
        x += Math.exp(-d2 / 60);
      }
      v[i * b + j] = x;
    }
  return v;
}

const letters = (n: number, seed: number) => {
  const r = new PCG(seed);
  let s = '';
  for (let i = 0; i < n; i++) s += 'ACGT'[r.int() % 4];
  return s;
};

/** Two chromosomes' sections (a block each), a grid following both: x chr1's, y chr4's (Juicebox's "chr1 vs chr4"). */
export function twoChromosomes(): Node {
  const b = createSeqSpec({ title: 'Contacts between two chromosomes' });
  const c1 = b.data.sequence('chr1', { alphabet: 'dna', letters: letters(ContactBins.a, 1) });
  const c4 = b.data.sequence('chr4', { alphabet: 'dna', letters: letters(ContactBins.b, 4) });
  const g = b.data.grid('contacts', { x: c1, y: c4, values: contactValues() });
  b.views.view('main', (v) => {
    v.section('chr1', {}, (s) => {
      s.block('chr1', { of: c1, label: 'chr1' });
      s.track('ruler').layer('ruler', { type: 'ruler' });
    });
    v.section('chr4', {}, (s) => {
      s.block('chr4', { of: c4, label: 'chr4' });
      s.track('ruler').layer('ruler', { type: 'ruler' });
    });
    v.section('contacts', { layout: 'grid', follow: '/views/main/chr1', follow_y: '/views/main/chr4' } as never)
      .track('contacts', { label: 'Contacts' })
      .layer('contacts', { type: 'grid', data: g } as never);
  });
  return b.build();
}

/** The kinship view's samples: three families, listed interleaved (a family's members apart in data order). */
export const Samples: readonly string[] = [
  'F2-a',
  'F1-a',
  'F3-a',
  'F1-b',
  'F2-b',
  'F3-b',
  'F1-c',
  'F3-c',
  'F2-c',
  'F1-d',
  'F2-d',
  'F3-d',
];
/** The genotyped sites. */
export const Sites = 160;

/**
 * The samples' genotypes (0, 1, 2 copies of the alternate allele), item-major (a site, then each sample): a family's
 * members share its founder's genotype at a site with probability 0.85, else draw their own. Seeded.
 */
export function genotypeCodes(): Uint8Array {
  const r = new PCG(29);
  const n = Samples.length;
  const out = new Uint8Array(Sites * n);
  for (let s = 0; s < Sites; s++) {
    const founders = [r.int() % 3, r.int() % 3, r.int() % 3];
    for (let k = 0; k < n; k++) {
      const family = Number(Samples[k][1]) - 1;
      out[s * n + k] = r.float() < 0.85 ? founders[family] : r.int() % 3;
    }
  }
  return out;
}

/** Identity by state over the sites: for samples `i`, `j` the mean of `1 − |gᵢ − gⱼ| / 2`, the grid's value. */
export function kinshipValues(codes: Uint8Array = genotypeCodes()): Float32Array {
  const n = Samples.length;
  const v = new Float32Array(n * n);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let s = 0; s < Sites; s++) sum += 1 - Math.abs(codes[s * n + i] - codes[s * n + j]) / 2;
      v[i * n + j] = sum / Sites;
    }
  return v;
}

/** The samples' rows system: the genotype matrix's states. */
export const SamplesSystem = '/data/gt/(states)';

/**
 * The genotype matrix (a row a sample, in data order) and a kinship grid over the same samples on both axes, ordered by
 * label (`member_order`): the families gather into blocks along the diagonal.
 */
export function kinship(): Node {
  const b = createSeqSpec({ title: 'Kinship beside the genotypes' });
  const chr = b.data.sequence('chr', { alphabet: 'dna', letters: letters(Sites, 9) });
  b.data.matrix('gt', { target: chr, codes: genotypeCodes(), states: [...Samples] } as never);
  const k = b.data.grid('kinship', { x: SamplesSystem, y: SamplesSystem, values: kinshipValues() } as never);
  b.views.view('main', (v) => {
    v.section('genotypes', { of: chr }, (s) => {
      s.track('ruler').layer('ruler', { type: 'ruler', units: 'bases' } as never);
      s.track('gt', { label: 'Genotypes', height: 160 }).layer('gt', { type: 'matrix', data: '/data/gt' } as never);
    });
    v.section('kinship', {
      layout: 'grid',
      of: SamplesSystem,
      of_y: SamplesSystem,
      member_order: { field: 'label' },
      member_order_y: { field: 'label' },
    } as never)
      .track('kinship', { label: 'Kinship (IBS)' })
      .layer('kinship', { type: 'grid', data: k } as never);
  });
  return b.build();
}
