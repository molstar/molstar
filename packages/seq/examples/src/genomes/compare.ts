/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import { createSeqSpec, genome, theme } from '@molstar/seqx-core/lang';
import { assemblyA, assemblyA2, assemblyB, type Files, fileSource, GenomesFiles, inlineA, inlineB } from './toys.js';

// The Genomes dataset's views of two genomes (examples.md, Genomes, The views): toy A over toy B with their alignment
// as ribbons (an inversion crossed, a translocation, a gap), the same alignment as a dot plot (two levels of grid
// tiles under the blocks), the contact map over toy A, and the liftover of A's genes onto A′ through the map between
// them. Each reads the committed files (`files`, the dataset's by default); the assemblies' lengths are their
// chrom.sizes files'.

/** The sections' addresses (what the tests read). */
export const ASection = '/views/main/toy-a';
export const BSection = '/views/main/toy-b';
export const A2Section = '/views/main/toy-a2';
export const DotsSection = '/views/main/dots';
export const ContactsSection = '/views/main/contacts';
export const SyntenyBridge = '/views/main/a-b';
export const LiftoverBridge = '/views/main/a-a2';

/**
 * Toy A above, toy B below, a bridge between them drawing the alignment's blocks as ribbons (sided pairs from the
 * `seq` file, read by region from either end, EC-208): A's chr1 inversion crossed (its blocks read reversed on B,
 * EC-206), A's chr3 tail carried to B's `2` (a translocation), chr1's run of N no ribbon (a gap).
 */
export function synteny(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Toy A over toy B' });
  const a = assemblyA(b, files);
  const g = assemblyB(b, files);
  const aln = b.data.pairs('alignment', { a, b: g, source: fileSource(b, 'alignment-file', 'seq', files.alignment) });
  b.views.view('main', (v) => {
    genome(v, 'toy-a', a);
    v.bridge('a-b', { a: ASection, b: BSection, label: 'A to B' } as never, (br) =>
      br.layer('blocks', { type: 'pairs', data: aln, legend: 'none' } as never),
    );
    genome(v, 'toy-b', g);
  });
  return b.build();
}

/**
 * The dot plot of toy A (across) against toy B (down), a grid section following both genome sections (EC-209): the
 * alignment's positions counted in cells of two levels (grid tiles, EC-223) and its blocks on top as diagonals (the
 * inversion an anti-diagonal in chr1 × 1, the translocation in chr3 × 2, chr1's run of N a break). A cell is two
 * positions (EC-145).
 */
export function dotPlot(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Toy A against toy B, a dot plot' });
  // inline: the grid follows both sections, whose blocks it needs before any file opens; names matched exactly (one
  // file lists both genomes' sequences, and B's `1` is A's `chr1` by the `chr` rule)
  const a = inlineA(b, 'toy-a', 'exact');
  const g = inlineB(b, 'toy-b', 'exact');
  const aln = b.data.pairs('alignment', { a, b: g, source: fileSource(b, 'alignment-file', 'seq', files.alignment) });
  const dots = b.data.grid('dots', { x: a, y: g, source: fileSource(b, 'dots-file', 'seq', files.dots) } as never);
  b.views.view('main', (v) => {
    genome(v, 'toy-a', a, { include: ['chr1', 'chr2', 'chr3'] });
    genome(v, 'toy-b', g, { include: ['1', '2', '3'] });
    v.section('dots', { layout: 'grid', follow: ASection, follow_y: BSection } as never).track(
      'dots',
      { label: 'Aligned positions' },
      (t) => {
        t.layer('cells', {
          type: 'grid',
          data: dots,
          color: theme.field('value', { palette: theme.continuous({ colors: 'blues' }) }),
        } as never);
        t.layer('blocks', { type: 'pairs', data: aln, legend: 'none' } as never);
      },
    );
  });
  return b.build();
}

/**
 * The Hi-C over toy A's chr1 to chr3 (EC-223): a grid section following the genome's, read by block pair from the
 * contact map's tiles (two levels; the upper triangle served, the lower its mirror), the chromosome gaps on both axes;
 * a hotspot between chr1 and chr3.
 */
export function contactMap(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Hi-C over toy A' });
  const a = inlineA(b);
  const contacts = b.data.grid('contacts', { x: a, source: fileSource(b, 'hic-file', 'seq', files.contacts) } as never);
  b.views.view('main', (v) => {
    genome(v, 'toy-a', a, { include: ['chr1', 'chr2', 'chr3'] });
    v.section('contacts', { layout: 'grid', follow: ASection } as never).track('contacts', {}, (t) =>
      t.layer('contacts', {
        type: 'grid',
        data: contacts,
        // counts on a log palette (a count of 1 to the peak's hundreds), an empty cell (0) no contact: left clear
        color: theme.field('value', { type: 'log', palette: theme.continuous({ colors: 'magma' }) }),
      } as never),
    );
  });
  return b.build();
}

/**
 * Toy A above with its genes, the map to toy A′ as a bridge's ribbons (the map's runs a `seq` source's sided pairs,
 * read by region: `map.source`, EC-211), A′ below with A's genes drawn through the map (EC-212): chr1's duplicated
 * 10 kb fanning to two places (held twice, side by side), chr2 shifted after its insertion, chr3's lost 4 kb no ribbon.
 */
export function liftover(files: Files = GenomesFiles): Node {
  const b = createSeqSpec({ title: 'Toy A lifted to A′' });
  const a = assemblyA(b, files);
  const a2 = assemblyA2(b, files);
  const chain = fileSource(b, 'liftover-file', 'seq', files.liftover);
  const map = b.data.map('a-a2', { from: a, to: a2, method: 'segments', source: chain } as never);
  const genes = b.data.ranges('genes', {
    target: a,
    source: fileSource(b, 'genes-file', 'gff3', files.gff3, files.gff3Index, { features: 'genes' }),
  } as never);
  const geneLayer = { type: 'ranges', data: genes, line: 'gene', legend: 'none' } as never;
  b.views.view('main', (v) => {
    genome(v, 'toy-a', a, { include: ['chr1', 'chr2', 'chr3'] }, (s) =>
      s.track('genes', { label: 'Genes (A)' }).layer('genes', geneLayer),
    );
    v.bridge('a-a2', { a: ASection, b: A2Section, label: 'A to A′' } as never, (br) =>
      br.layer('runs', { type: 'pairs', data: map, legend: 'none' } as never),
    );
    genome(v, 'toy-a2', a2, { include: ['chr1', 'chr2', 'chr3'] }, (s) =>
      s.track('genes', { label: 'Genes (lifted from A)' }).layer('genes', geneLayer),
    );
  });
  return b.build();
}
