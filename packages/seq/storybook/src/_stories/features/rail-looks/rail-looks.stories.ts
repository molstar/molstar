/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { railLooksStory } from './rail-looks-root.js';

const meta: Meta = { title: 'Design/Rail looks' };
export default meta;

/**
 * Three looks of the one rail (`navigator.rail.look`, internal: the story's alone until one is picked) beside today's, on
 * a ruler rail (zoomed in and at the fit), an alignment's rows rail and a lanes band (GFF3, TP53), in the light and the
 * dark theme: A a capsule, B the capsule in a faint groove, C the groove quiet until hovered. Each is live: hover and drag
 * the rail, its ends and the space beside the window.
 */
export const Looks: StoryObj = { render: () => railLooksStory() };
