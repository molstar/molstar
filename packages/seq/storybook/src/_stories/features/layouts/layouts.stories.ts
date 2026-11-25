/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { tracksPastTheCapStory } from './layouts-root.js';

const meta: Meta = { title: 'Features/Layouts' };
export default meta;

/**
 * A tracks section of a ruler, a lanes band and 30 genes, taller than its box: at the default cap (`section.maxPx`,
 * 600 px), a host's 300 and the section's own `height` of 200, light and dark, the canvas beside its SVG export. Scroll
 * it by the wheel, the rail at its right edge, or Page and Home/End with the pointer in it; the ruler stays.
 */
export const TracksPastTheCap: StoryObj = { render: () => tracksPastTheCapStory() };
