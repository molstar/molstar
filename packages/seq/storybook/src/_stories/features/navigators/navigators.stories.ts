/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

// Navigators as such (the navigators design): the ruler rail every linear section has by default, and a declared strip
// that replaces it; on p53, whose tracks show both best. The datasets' figures with navigators are their examples.

const meta: Meta = { title: 'Features/Navigators' };
export default meta;

const view = (id: string) => P53Views.find((v) => v.id === id)!;

/**
 * The ruler rail, every linear section's navigator by default: the whole axis on the ruler's bottom border, the window a
 * segment with a grip each end. Drag on the ruler to pan, drag a grip to zoom, click beside the window to centre there,
 * drag beside it to sweep a region, double-click to fit (fitted, to zoom back in). It grows while used, and is a quiet
 * hairline when the view shows everything at its zoom limit. An export leaves it out.
 */
export const Rail: StoryObj = { render: () => exampleStory(view('tracks')) };

/**
 * A declared strip over the same tracks replaces the rail: its picture the tracks small, its window the view's. Its
 * window is a slider: focus it and use the arrows, Shift+arrows, Home and End; its value is said when the view rests.
 * The SVG includes it, its window an outline titled with its positions.
 */
export const Strip: StoryObj = { render: () => exampleStory(view('tracks-strip'), { navigators: 'include' }) };
