/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec } from '@molstar/seqx-core/lang';
import { addDomains, addSecondary, addVariants, conservation, GeneExons, geneLetters, P53 } from './data.js';

// The small layers' views (the Features stories of the layers): the ruler, secondary structure, variants, values as a
// line and an area, ranges drawn as lines and with arrowheads, each on p53 and, where the axis matters, on its gene.

type Builder = ReturnType<typeof createSeqSpec>;

function data(b: Builder) {
  const p53 = b.data.sequence('p53', { alphabet: 'protein', letters: P53 });
  const gene = b.data.sequence('tp53', { alphabet: 'dna', letters: geneLetters() });
  b.data.map('p53-on-gene', { from: p53, to: gene, method: 'codon', exons: GeneExons, strand: '-' });
  return { p53, gene };
}

/** The ruler on the protein and, reversed on its coding pieces, on its gene: it counts down there. */
export function ruler(): Node {
  const b = createSeqSpec({ title: 'p53 rulers' });
  const { p53, gene } = data(b);
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    protein.track('letters', { label: 'p53' }).layer('letters', { type: 'letters' });
    const g = v.section('gene', { of: gene, direction: 'reversed', pieces: 'coding', follow: protein });
    g.track('ruler', { label: 'TP53' }).layer('ruler', { type: 'ruler' });
    g.track('letters', { label: 'TP53' }).layer('letters', { type: 'letters', complement: 'auto' });
  });
  return b.build();
}

/** Secondary structure on the protein and, through the codon map, on the gene that follows it. */
export function secondary(): Node {
  const b = createSeqSpec({ title: 'p53 secondary structure' });
  const { p53, gene } = data(b);
  const ss = addSecondary(b, p53);
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    protein.track('secondary', { label: 'Secondary' }).layer('secondary', { type: 'secondary', data: ss });
    const g = v.section('gene', { of: gene, follow: protein });
    g.track('secondary', { label: 'On the gene' }).layer('secondary', { type: 'secondary', data: ss });
  });
  return b.build();
}

/** The variants on the protein (two at R248, two at R273) and on its gene. */
export function variants(): Node {
  const b = createSeqSpec({ title: 'p53 variants' });
  const { p53, gene } = data(b);
  const vs = addVariants(b, p53);
  const domains = addDomains(b, p53);
  b.views.view('main', (v) => {
    const protein = v.section('protein', { of: p53 });
    protein.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    protein.track('variants', { label: 'Variants' }).layer('variants', { type: 'variants', data: vs });
    protein.track('domains', { label: 'Domains' }).layer('ranges', { type: 'ranges', data: domains });
    const g = v.section('gene', { of: gene, follow: protein });
    g.track('variants', { label: 'On the gene' }).layer('variants', { type: 'variants', data: vs });
  });
  return b.build();
}

/** The conservation as bars, a line and an area (two residues missing: the line breaks there). */
export function valuesLine(): Node {
  const b = createSeqSpec({ title: 'p53 conservation as a line and an area' });
  const { p53 } = data(b);
  const c = b.data.values('conservation', { target: p53, values: conservation() });
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    s.track('bars', { label: 'Bars' }).layer('bars', { type: 'values', data: c });
    s.track('line', { label: 'Line' }).layer('line', { type: 'values', data: c, encoding: 'line' });
    s.track('area', { label: 'Area' }).layer('area', { type: 'values', data: c, encoding: 'area' });
  });
  return b.build();
}

/** The domains drawn in each line style, with arrowheads on their ends where asked. */
export function rangesLines(): Node {
  const b = createSeqSpec({ title: 'p53 domains as lines and arrows' });
  const { p53 } = data(b);
  const domains = addDomains(b, p53);
  b.views.view('main', (v) => {
    const s = v.section('protein', { of: p53 });
    s.track('ruler', { label: 'Position' }).layer('ruler', { type: 'ruler' });
    for (const line of ['solid', 'dashed', 'dotted', 'wavy', 'double', 'box', 'highlight'])
      s.track(line, { label: line, height: 26 }).layer(line, { type: 'ranges', data: domains, line });
    s.track('arrows', { label: 'Arrows', height: 26 }).layer('arrows', { type: 'ranges', data: domains, arrow: 'end' });
    s.track('both', { label: 'Both', height: 26 }).layer('both', {
      type: 'ranges',
      data: domains,
      line: 'solid',
      arrow: 'both',
    });
  });
  return b.build();
}
