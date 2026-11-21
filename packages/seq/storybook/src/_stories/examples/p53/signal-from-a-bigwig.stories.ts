/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'signal-from-a-bigwig')!;

const meta: Meta = { title: 'Examples/p53/Signal from a bigWig' };
export default meta;

/**
 * The conservation signal over the gene's region from a bigWig (the opt-in `bigwig` runtime), read by region and zoom
 * level with a delay: pan and zoom to see coarse bins give way to a value a base.
 */
export const SignalFromABigWig: StoryObj = { render: () => exampleStory(view) };
