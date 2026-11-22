/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RegionReadsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { openHic } from './openings.js';

const view = (id: string) => RegionReadsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Region reads/Hi-C' };
export default meta;

/**
 * A generated Hi-C of one toy chromosome (384 positions, symmetric counts served at two levels), read in tiles by the
 * plot: the coarse tile arrives first and the fine ones after, a pending tile drawn from the coarser one held (blocky,
 * then sharp) with the plot’s box fixed. The plot opens on positions 100 to 250. The map is symmetric (the source says
 * so): a cell and its mirror read the same value. Hover a cell: its bin on each axis and its value (EC-169, EC-223).
 */
export const ByRegion: StoryObj = {
  render: () => exampleStory(view('hic-region'), { stacked: true, opening: openHic }),
};

/**
 * The Hi-C of a toy assembly’s six chromosomes as one system (their positions end to end, 1,376 bins), at the fit: the
 * source serves one level, the store’s summaries make the others and the whole-assembly level, so one read draws the
 * whole map. Zoom in and the finer levels are read by tiles. (Over the assembly’s blocks, with the chromosome gaps and
 * labels on both axes, waits for data over an assembly read by source.)
 */
export const OverTheGenome: StoryObj = { render: () => exampleStory(view('hic-genome'), { stacked: true }) };
