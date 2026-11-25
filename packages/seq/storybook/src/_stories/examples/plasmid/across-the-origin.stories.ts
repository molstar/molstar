/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlasmidViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const meta: Meta = { title: 'Examples/Plasmid/Features across the origin' };
export default meta;

/**
 * The rotated variant (its origin inside bla) with a bookmark on bla (the minus strand) and one on a generated primer
 * site (the plus strand), the two features that cross the origin. Click one: a single item in two runs, on the ring and
 * in the detail. Focusing a bookmark moves the view to its shortest arc, over the origin (EC-12); a Shift+drag selection
 * across the origin is pending.
 */
export const FeaturesAcrossTheOrigin: StoryObj = {
  render: () => exampleStory(PlasmidViews.find((v) => v.id === 'crossing')!),
};
