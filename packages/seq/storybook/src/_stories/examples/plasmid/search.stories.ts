/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlasmidViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';
import { AcrossTheOrigin } from './search.js';

const meta: Meta = { title: 'Examples/Plasmid/Search' };
export default meta;

/**
 * The ring with its linear detail and the stories' search box: twelve letters across the origin are searched at once,
 * one hit in two runs (the ring's arc over the origin, the detail's two pieces), its address one range from the end to
 * the start. Next reveals it by the shortest arc. Search its reverse complement to find it on the minus strand (both
 * strands by default; Forward finds nothing), or a restriction site (GAATTC, EcoRI: one hit, the same on both strands).
 */
export const SearchAcrossTheOrigin: StoryObj = {
  render: () => exampleStory(PlasmidViews.find((v) => v.id === 'crossing')!, { search: { initial: AcrossTheOrigin } }),
};
