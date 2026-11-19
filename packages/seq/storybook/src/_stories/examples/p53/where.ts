/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { storyId } from '../../story-id.js';

// Where each p53 view's story is (examples-and-stories): the dataset's figures are Examples; a story whose point is how
// something works is a Feature, shown on whichever dataset suits it. The dataset's index lists both.

/** A story: its title (under `Seq* (next)`) and its export. */
export interface Where {
  readonly title: string;
  readonly export: string;
}

/** The p53 views' figures: the stories under Examples. */
export const Figures: Readonly<Record<string, Where>> = {
  tracks: { title: 'Examples/p53/Tracks', export: 'Tracks' },
  'protein-on-gene': { title: 'Examples/p53/Protein on its gene', export: 'ProteinOnItsGene' },
  'gene-by-region': { title: 'Examples/p53/Gene by region', export: 'GeneByRegion' },
  chromosome: { title: 'Examples/p53/Chromosome', export: 'Chromosome' },
};

/** The p53 views whose point is a feature: the Features stories that show them. */
export const FeatureStories: Readonly<Record<string, Where>> = {
  text: { title: 'Features/Text', export: 'Text' },
  'text-centered': { title: 'Features/Text', export: 'Centered' },
  'text-under-tracks': { title: 'Features/Text', export: 'UnderTracks' },
  'one-axis': { title: 'Features/Coordinates', export: 'ProteinAndGeneOnOneAxis' },
  'genomic-follower': { title: 'Features/Coordinates', export: 'GenomicFollower' },
  'reversed-gene': { title: 'Features/Coordinates', export: 'ReversedGene' },
  'gene-region-ideogram': { title: 'Features/Navigators', export: 'GeneRegionOnChr17' },
  'tracks-strip': { title: 'Features/Navigators', export: 'Strip' },
  'protein-gene-strip': { title: 'Features/Navigators', export: 'ProteinStripWithItsGene' },
  themes: { title: 'Features/Themes', export: 'Themes' },
  story: { title: 'Features/Story', export: 'Story' },
};

/** The Features that act on p53, with the one line each: what to try. */
export const FeatureTries: readonly { readonly where: Where; readonly title: string; readonly summary: string }[] = [
  {
    where: { title: 'Features/Editing', export: 'Editing' },
    title: 'Editing',
    summary: 'The tracks edited: add and hide a track, set a theme, undo and redo, drag a height to preview it.',
  },
  {
    where: { title: 'Features/Addressing/Bookmarks', export: 'Bookmarks' },
    title: 'Bookmarks and links',
    summary: 'The gene’s region with bookmarks and a text whose links move the view, on transcripts that arrive later.',
  },
  {
    where: { title: 'Features/Selection', export: 'ProteinOnItsGene' },
    title: 'Selection',
    summary: 'Drag on the gene to select nucleotides; the protein shows them projected, in canvas and SVG.',
  },
];

/** The id of a story of `where`. */
export const idOf = (where: Where): string => storyId(where.title, where.export);
