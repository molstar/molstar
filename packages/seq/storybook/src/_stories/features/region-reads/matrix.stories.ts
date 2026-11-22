/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RegionReadsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { openDeletion } from './openings.js';

const view = (id: string) => RegionReadsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Region reads/Matrix' };
export default meta;

/**
 * 40 generated samples over a chromosome’s variants (a row a sample, a mark a variant at its own extent, coloured by
 * the genotype: none, one copy, both, uncalled in the missing colour), opened round one long deletion that six samples
 * carry: one solid bar across the SNVs it covers, the short ones drawn over it. Zoom out and the read is the coarse
 * level: a bin’s majority genotype per sample, its share as the colour’s strength. Hover a cell: the variant (the
 * shortest there: a SNP over the deletion), the sample and the genotype (EC-220, EC-222).
 */
export const SamplesByPositions: StoryObj = {
  render: () => exampleStory(view('matrix-samples'), { stacked: true, opening: openDeletion }),
};

/**
 * The same 40 samples on a ring: each sample a thin ring inside the chromosome’s arc, a run at most a device column wide
 * a ring (the polar warp draws the image as runs, so the SVG is rects, not a picture). Wheel over the ring for the
 * window (EC-108, EC-220).
 */
export const OnARing: StoryObj = { render: () => exampleStory(view('matrix-ring'), { stacked: true }) };

/**
 * A methylation fraction per CpG and sample (eight samples, an island at 240 to 262 kb where it drops), about a tenth of
 * the values missing: a missing value is left clear, never zero. Zoomed out, a bin’s mean of the values it has (EC-221).
 */
export const MethylationWithMissingValues: StoryObj = {
  render: () => exampleStory(view('matrix-methylation'), { stacked: true }),
};
