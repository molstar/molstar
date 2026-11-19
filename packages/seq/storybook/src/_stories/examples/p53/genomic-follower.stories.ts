/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'genomic-follower')!;

const meta: Meta = { title: 'Examples/p53/Genomic follower' };
export default meta;

/**
 * For contrast: the gene follows its protein in genomic order, its introns shown. A pan of the protein moves the gene
 * the other way (minus strand); in Protein on its gene both move together.
 */
export const GenomicFollower: StoryObj = { render: () => exampleStory(view) };
