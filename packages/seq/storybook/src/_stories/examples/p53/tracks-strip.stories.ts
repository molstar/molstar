/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'tracks-strip')!;

const meta: Meta = { title: 'Examples/p53/Tracks with a strip' };
export default meta;

/**
 * p53's tracks with a strip of them small above (the domains in their lanes): drag the window to pan, its grips to zoom,
 * click outside it to centre there, drag outside it to sweep a region; Alt+Left goes back. The ruler has no rail: the
 * strip replaces it.
 */
export const TracksWithAStrip: StoryObj = { render: () => exampleStory(view, { navigators: 'include' }) };
