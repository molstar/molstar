/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import * as addressing from './p53/addressing.js';
import * as axis from './p53/axis.js';
import * as bookmarks from './p53/bookmarks.js';
import * as geneByRegion from './p53/gene-by-region.js';
import * as proteinOnGene from './p53/protein-on-gene.js';
import * as tracks from './p53/tracks.js';

// The examples' views by dataset: what a dataset's index lists and a test walks.

export interface ExampleView {
  readonly dataset: string;
  /** The module's name under the dataset (`tracks`). */
  readonly id: string;
  readonly title: string;
  /** What the view shows, a sentence. */
  readonly summary: string;
  build(): Node;
  /**
   * What a viewer of the view needs besides the default registry (a runtime of a format its sources name); answers
   * delayed by `delayMs` (the view's own pacing when absent; 0 in tests: no timer).
   */
  registry?(options?: { readonly delayMs?: number }): readonly SeqRegistryEntry[];
}

export const P53Views: readonly ExampleView[] = [
  {
    dataset: 'p53',
    id: 'tracks',
    title: 'Tracks',
    summary: 'The protein’s letters, its domains in lanes and a conservation signal with two residues missing.',
    build: tracks.build,
  },
  {
    dataset: 'p53',
    id: 'protein-on-gene',
    title: 'Protein on its gene',
    summary:
      'The gene’s section follows the protein’s through the codon map (minus strand, three exons), reading along it (`direction: follow`) on its coding pieces: each intron a hairline, three nucleotides a residue, both moving the same way.',
    build: proteinOnGene.build,
  },
  {
    dataset: 'p53',
    id: 'one-axis',
    title: 'Protein and gene on one axis',
    summary:
      'One gene section with the protein’s letters a codon wide, its domains and the gene’s letters: it reads along the protein (`{ along: p53 }`) on its coding pieces; a codon split by an intron is drawn in each part; R175H, the gene’s one variant, is marked as a conflict.',
    build: axis.oneAxis,
  },
  {
    dataset: 'p53',
    id: 'genomic-follower',
    title: 'Genomic follower',
    summary:
      'For contrast: the gene follows its protein in genomic order (`direction: genomic`, introns shown), so on the minus strand it moves against the protein.',
    build: axis.genomicFollower,
  },
  {
    dataset: 'p53',
    id: 'reversed-gene',
    title: 'Reversed gene',
    summary:
      'The gene alone, reversed: positions count down left to right, its letters read on the other strand (the complement), the domains through the codon map; Shift+C collapses its introns.',
    build: axis.reversedGene,
  },
  {
    dataset: 'p53',
    id: 'gene-by-region',
    title: 'Gene by region',
    summary:
      'The gene’s region read region by region: transcripts arrive as dynamic nodes, a track each, and TP53’s bring their domains, shown through each transcript’s codon map.',
    build: geneByRegion.build,
    registry: geneByRegion.registry,
  },
];

/** Views made for a feature's stories (`Features/…`), over the same datasets. */
export const FeatureViews: readonly ExampleView[] = [
  {
    dataset: 'p53',
    id: 'addressing',
    title: 'Addressing',
    summary: 'The tracks view with refs naming its domains and hotspots, and a text whose links act on addresses.',
    build: addressing.build,
  },
  {
    dataset: 'p53',
    id: 'bookmarks',
    title: 'Bookmarks',
    summary: 'Gene by region with bookmarks: refs on a range and on transcripts that arrive with their region.',
    build: bookmarks.build,
    registry: bookmarks.registry,
  },
];

export const ExampleViews: readonly ExampleView[] = [...P53Views, ...FeatureViews];
