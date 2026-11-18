/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { selectionStory } from './selection-root.js';

// Selection: the default grammar on the two p53 views (the interaction grammar); the story's root and what it does are
// in `selection-root.ts`.

const meta: Meta = { title: 'Features/Selection' };
export default meta;

const story = (id: string) => selectionStory(P53Views.find((v) => v.id === id)!);

/** p53's tracks: select residues, domains and ranges; the cursor; the SVG beside it draws the same selection. */
export const Tracks: StoryObj = { render: () => story('tracks') };

/** The protein on its gene: what is selected or hovered on one section shows on the other, lighter. */
export const ProteinOnItsGene: StoryObj = { render: () => story('protein-on-gene') };
