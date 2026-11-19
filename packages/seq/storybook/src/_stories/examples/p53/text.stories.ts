/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => P53Views.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/p53/Text' };
export default meta;

/**
 * The protein as text: lines of whole groups of ten wrapped to the width, a number a line and a group, its secondary
 * structure and domains as thin bands under the letters (a band left out of a line it is empty on). Resize: the lines
 * reflow with the first line kept at the top.
 */
export const Text: StoryObj = { render: () => exampleStory(view('text')) };

/** Centered in at most 640 px, the group numbers inline in the gaps, a number at each end of a line. */
export const Centered: StoryObj = { render: () => exampleStory(view('text-centered')) };

/** A tracks overview of the domains and the protein as text under it: two layouts in one view. */
export const UnderTracks: StoryObj = { render: () => exampleStory(view('text-under-tracks')) };
