/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'tracks')!;

const meta: Meta = { title: 'Examples/p53/Tracks' };
export default meta;

/** p53's letters, its domains in lanes and a conservation signal: scroll and zoom the canvas, hover a domain. */
export const Tracks: StoryObj = { render: () => exampleStory(view) };
