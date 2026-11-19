/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'protein-gene-strip')!;

const meta: Meta = { title: 'Examples/p53/Protein strip with its gene' };
export default meta;

/**
 * The protein on its gene with a strip of the protein: its window and its gene follower's, carried along the codon map,
 * drawn together; a drag on the gene's window moves the protein, which leads.
 */
export const ProteinStripWithItsGene: StoryObj = { render: () => exampleStory(view, { navigators: 'include' }) };
