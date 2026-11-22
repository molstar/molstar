/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AcrossViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { openZoomed } from './openings.js';

const view = (id: string) => AcrossViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Across sections' };
export default meta;

/**
 * Peach above, grape below, a 96 px band between (generated genomes; the alignment’s blocks as pairs in the `seq`
 * format, read as the sides need them): each block a ribbon filled at 0.35 in its peach chromosome’s hue; Pp02’s
 * inversion crossed, Pp03’s tail carried to grape’s chr4 (a translocation), Pp01’s lost stretch no ribbon (a gap).
 * Drag or wheel over the band: both genomes move; Shift+A at a side brings the other to what it shows (EC-206, EC-218).
 */
export const TwoGenomesAtTheFit: StoryObj = { render: () => exampleStory(view('across-fit'), { stacked: true }) };

/**
 * The same genomes aligned in about 6,000 blocks: past the 4,096 a band draws, bin pairs of each side’s width, their
 * opacity by the bases aligned; zoom in past the gate and the blocks replace them with no gap (EC-129).
 */
export const ZoomedOutBins: StoryObj = { render: () => exampleStory(view('across-bins'), { stacked: true }) };

/**
 * Peach zoomed into Pp02 from 1.2 to 1.5 Mb, grape at the fit: the zoomed side’s blocks fan to their true x on grape,
 * narrowing to grape’s 300 kb; ends past a side’s edge leave the band sideways.
 */
export const OneSideZoomedIn: StoryObj = {
  render: () => exampleStory(view('across-zoomed'), { stacked: true, opening: openZoomed }),
};

/**
 * Grape without its chr4 (its `include`): the ribbons of peach Pp03’s tail have no place below, each a quarter-height
 * taper from its peach end, never an issue.
 */
export const AStub: StoryObj = { render: () => exampleStory(view('across-stub'), { stacked: true }) };

/**
 * One chromosome each, grape’s reads answered after peach’s: the pending hatch along the edge of the side not loaded,
 * then the ribbons (EC-209).
 */
export const PendingASideAtATime: StoryObj = {
  render: () => exampleStory(view('across-pending'), { stacked: true }),
};

/** The bridge hidden (kept in the spec): the two genomes stacked with the sections’ gap, the band gone. */
export const AReaderWithoutBridges: StoryObj = {
  render: () => exampleStory(view('across-hidden'), { stacked: true }),
};

/**
 * Peach above and below, its duplications drawn both ways (`mirror`): one within Pp02, one onto Pp03, one reversed
 * onto Pp04; the lower section follows the upper, so the two move together.
 */
export const OneGenomeAgainstItself: StoryObj = {
  render: () => exampleStory(view('across-self'), { stacked: true }),
};

/**
 * Peach, grape and plum, a bridge between each two: grape’s chr1 tail carried reversed onto plum’s Pd4, its chr3 and
 * chr4 fused in plum’s Pd3. Shift+A on grape brings both of its neighbours over in one step (EC-235).
 */
export const ThreeGenomesChained: StoryObj = {
  render: () => exampleStory(view('across-chain'), { stacked: true }),
};
