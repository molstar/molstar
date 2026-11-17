/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'protein-on-gene')!;

const meta: Meta = { title: 'Examples/p53/Protein on its gene' };
export default meta;

/** The gene's section follows the protein's through the codon map; its letters are read from a fixture file. */
export const ProteinOnItsGene: StoryObj = { render: () => exampleStory(view) };
