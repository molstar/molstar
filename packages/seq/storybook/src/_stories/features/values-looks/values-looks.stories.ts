/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { valuesLooksStory } from './values-looks-root.js';

const meta: Meta = { title: 'Design/Values looks' };
export default meta;

/**
 * The values layer's look: today's beside three variants (a stepped profile where a value is under 4 px, bars from it;
 * separate rounded bars; a quieter fill) on p53's conservation and a bigWig signal, zoomed out and in, in the light and
 * the dark theme. Each is live: hover a value to light it.
 */
export const Looks: StoryObj = { render: () => valuesLooksStory() };
