/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { AcrossViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { openZoomed, openZoomedGrape } from './openings.js';

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
 * narrowing to grape’s 300 kb; grape reads its own blocks, so every block is drawn once, those anchored in peach’s
 * window by peach, the rest by grape, leaving the band sideways.
 */
export const OneSideZoomedIn: StoryObj = {
  render: () => exampleStory(view('across-zoomed'), { stacked: true, opening: openZoomed }),
};

/**
 * The other way, over the dense alignment: grape zoomed into chr2 from 1.2 to 1.5 Mb, its blocks; peach at the fit,
 * its bin pairs. Each block is drawn once: grape’s window (widened to peach’s bins) by grape, the rest by peach’s bins.
 */
export const OneSideZoomedInTheOtherWay: StoryObj = {
  render: () => exampleStory(view('across-bins'), { stacked: true, opening: openZoomedGrape }),
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

/**
 * Peach across, grape down: a grid over two sections of blocks, a cell a pair of chromosomes (EC-209). Each alignment
 * block a diagonal; Pp02’s inversion an anti-diagonal in Pp02 × chr2, Pp03’s tail in Pp03 × chr4 (a translocation),
 * Pp01’s lost stretch a break in its diagonal. Pan either genome and its axis of the plot follows.
 */
export const TwoGenomesAsADotPlot: StoryObj = {
  render: () => exampleStory(view('across-dot-plot'), { stacked: true }),
};

/**
 * chr1’s section across, chr4’s down (a block each), their contacts a grid following both (EC-223): two hotspots where
 * the chromosomes touch, over a low floor.
 */
export const ContactsBetweenTwoChromosomes: StoryObj = {
  render: () => exampleStory(view('across-contacts'), { stacked: true }),
};

/**
 * A liftover (EC-212): Rose v1 above with its genes, the chain to v2 as ribbons (a map read by region, its runs a sided
 * pairs source’s rows), v2 below with v1’s genes drawn through the same chain. chr1’s 1.2–2.0 Mb fans to two places (a
 * duplication: its gene drawn twice, EC-20), chr2’s lost stretch has no ribbon and its gene is not drawn (EC-21), chr2’s
 * inversion crosses, chr3’s tail lands at the end of v2’s chr1. Hover a gene: its images light below, and the card of
 * the duplicated one lists its two places.
 */
export const Liftover: StoryObj = { render: () => exampleStory(view('across-liftover'), { stacked: true }) };

/**
 * Twelve samples’ genotypes as a matrix (a row a sample, in data order) and their kinship on the same samples as a
 * grid on categorical axes, ordered by label (`member_order`, EC-35): the three families gather into blocks along the
 * diagonal. The members’ labels run down the side and at 45° along the top; a click on one selects the sample.
 */
export const KinshipBesideTheGenotypes: StoryObj = {
  render: () => exampleStory(view('across-kinship'), { stacked: true }),
};
