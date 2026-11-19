/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LayerViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

// The small layers, a story each, on p53: the canvas and the SVG export of one spec side by side, in the light and the
// dark theme (the toolbar's theme switch redraws both with the other theme's tokens). Look at what the summary says.

const meta: Meta = { title: 'Features/Layers' };
export default meta;

const story = (id: string) => {
  const view = LayerViews.find((v) => v.id === id)!;
  return { render: () => exampleStory(view) } satisfies StoryObj;
};

/**
 * The ruler: labels at steps of 1, 2, 5 × 10ⁿ, none overlapping; on the gene (reversed, its coding pieces) they count
 * down. Zoom in and out: the step changes, a label that does not fit is skipped, never shrunk.
 */
export const Ruler: StoryObj = story('ruler');

/**
 * Secondary structure: helices as blocks, strands with an arrowhead on their C-terminal end, turns thin. On the gene
 * through the codon map the strands still point from N to C on the screen. Pan a strand's end out of view: its
 * arrowhead goes, never a stub.
 */
export const SecondaryStructure: StoryObj = story('secondary');

/**
 * Variants: a stem and a head a site; R248Q and R248W, R273H and R273C stand stacked on one stem. Zoom out: the heads
 * become a count strip in the density colors.
 */
export const Variants: StoryObj = story('variants');

/**
 * Values as bars, a line and an area: the line is held flat over half a cell at its ends and where residues 101 and
 * 102 are missing, where it breaks; the area is the line's, under it.
 */
export const ValuesAsALineAndAnArea: StoryObj = story('values-line');

/**
 * Ranges as lines: solid, dashed, dotted, wavy, double, box and highlight, and arrowheads on the domains' ends. Pan a
 * dashed domain: the dashes travel with it.
 */
export const RangesAsLinesAndArrows: StoryObj = story('ranges-lines');

/**
 * Tracks in their wells: each well a squircle (a third of its height, capped), its data clipped to it, so the letters'
 * cells lose a sliver at the corners. The letters are padded on top, the domains top and bottom: the data sits inside
 * the padding, and a hover or a hit there finds nothing. Hover a cell: its outline is a squircle too.
 */
export const TracksInTheirWells: StoryObj = story('wells');
