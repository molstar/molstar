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
 * The values layer's look on p53's conservation and a bigWig signal, zoomed out (a stepped profile) and zoomed in (bars),
 * in the light and the dark theme. Each is live: hover a value to light it.
 */
export const Looks: StoryObj = { render: () => valuesLooksStory() };
