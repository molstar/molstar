/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews, GlobinsViews, P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { openStrip } from './compose.js';

// Navigators as such (the navigators design): the ruler rail every linear section has by default, and a declared strip
// that replaces it; on p53, whose tracks show both best; a strip with a follower, and the ideogram of a chromosome over
// a gene's region; the vertical forms over the Globins' rows (the rows rail, the vertical overview). The other datasets'
// strips are drawn in their own examples.

const meta: Meta = { title: 'Features/Navigators' };
export default meta;

const view = (id: string) => P53Views.find((v) => v.id === id)!;

/**
 * The ruler rail, every linear section's navigator by default: the whole axis on the ruler's bottom border, the window a
 * segment with a grip each end. Drag on the ruler to pan, drag a grip to zoom, click beside the window to centre there,
 * drag beside it to sweep a region, double-click to fit (fitted, to zoom back in). It grows while used, and is a quiet
 * hairline when the view shows everything at its zoom limit. An export leaves it out.
 */
export const Rail: StoryObj = { render: () => exampleStory(view('tracks')) };

/**
 * A declared strip over the same tracks replaces the rail: its picture the tracks small, its window the view's. Its
 * window is a slider: focus it and use the arrows, Shift+arrows, Home and End; its value is said when the view rests.
 * The SVG includes it, its window an outline titled with its positions. It opens zoomed in on the DNA-binding domain, so
 * the window is smaller than the strip and can be dragged.
 */
export const Strip: StoryObj = {
  render: () =>
    exampleStory(view('tracks-strip'), {
      navigators: 'include',
      opening: openStrip,
    }),
};

/**
 * The protein on its gene with a strip of the protein: its window and its gene follower's, carried along the codon map,
 * drawn together; a drag on the gene's window moves the protein, which leads.
 */
export const ProteinStripWithItsGene: StoryObj = {
  render: () => exampleStory(view('protein-gene-strip'), { navigators: 'include' }),
};

/**
 * The gene's region with chromosome 17's ideogram above it (GRCh38's bands, the centromere two triangles): the region,
 * a sliver of the chromosome, marked on it, and a zoom strip under it over the region alone, joined to the mark by two
 * lines. Hover a band for its name, the mark for the region; drag the zoom strip's window, its grips, click outside it,
 * double-click. The SVG beside it includes the ideogram and its zoom strip, the window an outline titled with its
 * positions.
 */
export const GeneRegionOnChr17: StoryObj = {
  render: () => exampleStory(view('gene-region-ideogram'), { navigators: 'include' }),
};

/**
 * The rows rail, an alignment's rows' navigator when none is declared: the rail upright at the right edge of the rows'
 * band, beside their numbers, its window the rows in view. Drag it to scroll the rows, a grip to zoom them, click
 * beside the window to centre there; `]` and `[` zoom them. Idle while every row shows.
 */
export const RowsRail: StoryObj = {
  render: () => exampleStory(FeatureViews.find((v) => v.id === 'rows-rail')!),
};

/**
 * The vertical overview of the same rows, declared (it replaces the rows rail): every row summarized down a column at
 * the right, its window the rows and the columns in view, a grip on each of its four sides. Drag the window to move
 * both, its top and bottom grips to zoom the rows, its left and right grips to zoom the columns; click outside it to
 * centre there. The SVG includes it.
 */
export const VerticalOverview: StoryObj = {
  render: () => exampleStory(GlobinsViews.find((v) => v.id === 'family')!, { navigators: 'include' }),
};
