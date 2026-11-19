/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Node } from '@molstar/core/repr/tree';
import type { SeqRegistryEntry } from '@molstar/seqx-core/context';
import * as complexBlocks from './complex/blocks.js';
import * as plasmidRing from './plasmid/ring.js';
import * as addressing from './p53/addressing.js';
import * as axis from './p53/axis.js';
import * as bookmarks from './p53/bookmarks.js';
import * as geneByRegion from './p53/gene-by-region.js';
import * as layers from './p53/layers.js';
import * as proteinOnGene from './p53/protein-on-gene.js';
import * as story from './p53/story.js';
import * as themes from './p53/themes.js';
import * as text from './p53/text.js';
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
    id: 'text',
    title: 'Text',
    summary:
      'The protein as text: lines of whole groups of ten wrapped to the width, a number a line and a group, its secondary structure and domains as thin bands under the letters, each left out of a line it is empty on.',
    build: text.build,
  },
  {
    dataset: 'p53',
    id: 'text-centered',
    title: 'Text, centered',
    summary:
      'The text figure centered in at most 640 px, its group numbers inline in the gaps, a number at each end of a line.',
    build: text.centered,
  },
  {
    dataset: 'p53',
    id: 'text-under-tracks',
    title: 'Text under its tracks',
    summary:
      'A tracks overview of the domains and the protein as text under it, in one view: the text section’s own layout.',
    build: text.underTracks,
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
  {
    dataset: 'p53',
    id: 'themes',
    title: 'Themes',
    summary:
      'The protein’s letters under each scheme that fits it; domains by category with their legend; a confidence as pLDDT’s four bands and as a ramp; two signals on one shared scale; variants by consequence, sized by a log-scaled frequency; secondary structure. Its export draws the legends below.',
    build: themes.build,
  },
  {
    dataset: 'p53',
    id: 'story',
    title: 'A story',
    summary:
      'The tracks as a story of four steps: the whole protein, its DNA-binding domain (the conservation taller), the R175H hotspot (marked by a track of its own, the domains recolored) and the tetramerization domain (the conservation gone); each a view and a patch over the spec.',
    build: story.build,
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
  {
    dataset: 'p53',
    id: 'themes-editing',
    title: 'Editing themes',
    summary:
      'The themes view edited: a scheme previewed while hovered and picked in one patch; the shared scale’s domain set for both signals in one step and undone in one.',
    build: themes.build,
  },
];

/** Views made for the small layers' stories (`Features/Layers`), over p53. */
export const LayerViews: readonly ExampleView[] = [
  {
    dataset: 'p53',
    id: 'ruler',
    title: 'Ruler',
    summary:
      'The ruler on the protein and, reversed on its coding pieces, on its gene: the labels are steps of 1, 2 or 5 × 10ⁿ that never overlap, and the gene’s count down.',
    build: layers.ruler,
  },
  {
    dataset: 'p53',
    id: 'secondary',
    title: 'Secondary structure',
    summary:
      'Helices and strands as blocks, a turn thin; each strand ends in an arrowhead on the piece holding its C-terminal end, on the protein and through the codon map on the gene reading along it.',
    build: layers.secondary,
  },
  {
    dataset: 'p53',
    id: 'variants',
    title: 'Variants',
    summary:
      'A stem and a head a site, the heads of one site stacked (R248Q and R248W, R273H and R273C), on the protein and on the gene; zoom out and they become a count strip.',
    build: layers.variants,
  },
  {
    dataset: 'p53',
    id: 'values-line',
    title: 'Values as a line and an area',
    summary:
      'The conservation as bars, a line (held over half cells at its ends and at the two missing residues, broken there) and an area under it.',
    build: layers.valuesLine,
  },
  {
    dataset: 'p53',
    id: 'ranges-lines',
    title: 'Ranges as lines and arrows',
    summary:
      'The domains in each line style (the dash and the wave run from each domain’s start, so they move with it as you pan) and with arrowheads on their ends.',
    build: layers.rangesLines,
  },
];

/** The Plasmid (a stand-in until its dataset is recorded): the ring with its linear detail, the rotated linear view. */
export const PlasmidViews: readonly ExampleView[] = [
  {
    dataset: 'plasmid',
    id: 'ring',
    title: 'Ring with a linear detail',
    summary:
      'The construct as rings (each strand’s features, the sites, a ruler; labels outside, the name in the hole) and a linear detail following the ring’s window; bla and a primer site cross the origin, each one arc. Drag the window, or pan the detail: the other follows.',
    build: () => plasmidRing.build(),
  },
  {
    dataset: 'plasmid',
    id: 'rotated-linear',
    title: 'Rotated linear view',
    summary:
      'The construct in tracks from inside bla over the origin to past lacZα: a `range` past the origin, two pieces meeting at the seam with no break.',
    build: plasmidRing.rotatedLinear,
  },
];

/** The Complex (a stand-in until its dataset is recorded): its chains and ligands as blocks, three ways. */
export const ComplexViews: readonly ExampleView[] = [
  {
    dataset: 'complex',
    id: 'blocks',
    title: 'Chains as blocks',
    summary:
      'CDK2 and cyclin A in a row of blocks, polymers first, then ATP as a one-position block; the two copies of the complex merged, each block labeled ×2.',
    build: () => complexBlocks.build(),
  },
  {
    dataset: 'complex',
    id: 'blocks-split',
    title: 'Copies split',
    summary: 'The same blocks with the copies split (`copies: each`): a block a chain and a ligand.',
    build: complexBlocks.split,
  },
  {
    dataset: 'complex',
    id: 'blocks-stacked',
    title: 'Chains stacked',
    summary: 'A block a row of the tracks, aligned at their first residue, one zoom (`arrange: stack`).',
    build: complexBlocks.stacked,
  },
  {
    dataset: 'complex',
    id: 'blocks-ring',
    title: 'Chains on a ring',
    summary: 'The blocks on a ring, a block an arc with a gap between, their labels outside with the features’.',
    build: complexBlocks.ring,
  },
];

export const ExampleViews: readonly ExampleView[] = [
  ...P53Views,
  ...FeatureViews,
  ...LayerViews,
  ...PlasmidViews,
  ...ComplexViews,
];
