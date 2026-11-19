/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';
import { openOnRegion } from './opening.js';

const view = P53Views.find((v) => v.id === 'gene-by-region')!;

const meta: Meta = { title: 'Examples/p53/Gene by region' };
export default meta;

/**
 * The region read region by region with a delay: transcripts arrive as dynamic nodes, TP53's with their domains. It
 * opens on a region of the gene, so the toolbar's Back (to the whole gene) and Forward move; each is disabled, with its
 * reason as its tooltip, while there is no view to go to.
 */
export const GeneByRegion: StoryObj = { render: () => exampleStory(view, { opening: openOnRegion }) };
