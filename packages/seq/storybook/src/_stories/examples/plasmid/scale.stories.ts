/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlasmidViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const meta: Meta = { title: 'Examples/Plasmid/Ring at scale' };
export default meta;

/**
 * A generated 200 kb construct with 2,000 features on both strands, a tenth of them across the origin, as a ring with
 * its linear detail. The features show in the detail; on the ring, past its lanes (the 200 across the origin lie over
 * one position), each strand's features are its density, a band round the ring.
 */
export const RingAtScale: StoryObj = { render: () => exampleStory(PlasmidViews.find((v) => v.id === 'scale')!) };
