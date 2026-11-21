/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { GenomeViews } from './views.js';

// Genome coordinates (Design › Genome coordinates): chromosomes as blocks of one section, a direction and a scale per
// block, gaps by a share and small blocks in runs, two assemblies as groups. Toy assemblies generated in the story;
// each story the canvas and the SVG of one view, in the theme the toolbar picks (light or dark).

const meta: Meta = { title: 'Features/Genome Coordinates' };
export default meta;

const view = (id: string) => GenomeViews.find((v) => v.id === id)!;

/**
 * One assembly in a row: its 25 sequences end to end, each block its own ruler in Mb and its label at its start; the
 * gaps share at most a quarter of the width, 12 px each at most (chrM, under a pixel at the fit, appears as you zoom in).
 * Wheel to zoom into a chromosome: the pan carries across the gap into the next.
 */
export const OneAssemblyInARow: StoryObj = { render: () => exampleStory(view('row'), { stacked: true }) };

/**
 * The assembly as a ring, arcs clockwise from the top, labels outside, its name and length in the hole; the linear
 * detail under it has the same blocks, so it follows the ring's window (drag the rim's window, or pan the detail).
 */
export const OneAssemblyOnARing: StoryObj = { render: () => exampleStory(view('ring'), { stacked: true }) };

/** Two assemblies on one ring (peach and grape): a group each, a wider gap between them, an arc and a name outside. */
export const TwoAssembliesOnOneRing: StoryObj = {
  render: () => exampleStory(view('two-on-a-ring'), { stacked: true }),
};

/**
 * B's chromosomes in A's order: chr2 cut at an inversion's edges, its middle block reversed (`chr2 inv ←`, its ruler
 * counting down); chr4 scaled four times (`chr4 ×4`, its ruler at its own density).
 */
export const BInAsOrder: StoryObj = {
  name: 'B In A’s Order',
  render: () => exampleStory(view('b-in-a-order'), { stacked: true }),
};

/**
 * A draft of 24 chromosomes and 2,976 contigs, on a ring and in a row: the contigs narrower than 2 px are one run (a
 * band, its count where it has room); zoom into the run and its contigs part again.
 */
export const ThreeThousandContigs: StoryObj = {
  name: '3,000 Contigs',
  render: () => exampleStory(view('contigs'), { stacked: true }),
};
