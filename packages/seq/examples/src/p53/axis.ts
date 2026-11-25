/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, GeneExons, GeneFile, P53 } from './data.js';

// The axis's p53 views (the application drives the composition: each spells its direction and pieces): the protein
// and its gene on one axis, the genomic follower for contrast, the reversed gene with its complement. The gene's
// letters are read from the committed file, as the protein-on-gene view reads them.

type Builder = ReturnType<typeof createSeqSpec>;

function data(b: Builder) {
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const file = b.data.source('tp53-gene', { format: 'seq', uri: GeneFile });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', source: file });
  const map = b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GeneExons, strand: '-' });
  const domains = addDomains(b, p53);
  return { p53, gene, map, domains };
}

/**
 * The protein and its gene on one axis: a gene section with the protein's letters a codon wide (through the codon
 * map), its domains and the gene's letters; it reads along the protein (`{ along: p53 }`: a minus-strand gene laid out
 * N to C) and shows its coding pieces, each intron a hairline.
 */
export function oneAxis(): Node {
  const b = createSeqSpec({ title: 'p53 and its gene on one axis' });
  const { p53, gene, domains } = data(b);
  b.views.view('main', (v) => {
    const g = v.section('gene', { of: gene, direction: { along: p53 }, pieces: 'coding' });
    g.track('protein', { label: 'p53' }).layer('letters', { type: 'letters', data: p53 });
    g.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters' });
  });
  return b.build();
}

/**
 * The genomic follower, for contrast: the gene follows its protein in genomic order (`direction: genomic`, its introns
 * shown), so on the minus strand it moves against the protein.
 */
export function genomicFollower(): Node {
  const b = createSeqSpec({ title: 'p53 and its gene in genomic order' });
  const { p53, gene, domains } = data(b);
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    const g = v.section('gene', { of: gene, follow: protein, direction: 'genomic', pieces: 'whole' });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters' });
    g.track('domains', { label: 'Domains on the gene' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}

/**
 * The reversed gene: the gene alone, `reversed`, the ruler counting down; its letters read on the other strand (the
 * complement, `complement: auto`), the protein's domains through the map, the exons collapsed when `section.collapse`
 * asks (`collapsed: coding` is the default).
 */
export function reversedGene(): Node {
  const b = createSeqSpec({ title: 'TP53 reversed' });
  const { gene, domains } = data(b);
  b.views.view('main', (v) => {
    const g = v.section('gene', { of: gene, direction: 'reversed' });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters', complement: 'auto' });
    g.track('domains', { label: 'Domains on the gene' }).layer('ranges', { type: 'ranges', data: domains });
  });
  return b.build();
}
