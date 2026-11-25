/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlasmidViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const meta: Meta = { title: 'Examples/Plasmid/Rotated linear view' };
export default meta;

/**
 * The published record in tracks from inside bla over the origin to past lacZα: a `range` that starts after it ends,
 * two pieces meeting at the seam with no break, every feature whole. The range's start is the one param to set: the
 * Editing story does it as a patch.
 */
export const RotatedLinearView: StoryObj = {
  render: () => exampleStory(PlasmidViews.find((v) => v.id === 'rotated-linear')!),
};
