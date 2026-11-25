/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'chromosome')!;

const meta: Meta = { title: 'Examples/p53/Chromosome' };
export default meta;

/**
 * Chromosome 17 at its GRCh38 length (83 Mb), generated region by region and never stored, opened fit to the whole: the
 * genes a heat strip at first, then bars and items as you zoom; each gene's exons a track of its own once its region
 * arrives; a GC-like signal by levels; the letters when they fit; chr17's ideogram above with the window. Double-click
 * TP53 (7.67 Mb, minus strand) to focus it, zoom on to its letters (the p53 gene's); 0 fits again.
 */
export const Chromosome: StoryObj = { render: () => exampleStory(view, { navigators: 'include' }) };
