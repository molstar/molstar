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

/**
 * The gene's section follows the protein's through the codon map, reading along it on its coding pieces (each intron a
 * hairline, three nucleotides a residue); its letters are read from a fixture file. A pan, a zoom or a pinch on either
 * section moves both the same way: on the gene it moves the protein through the map read back.
 */
export const ProteinOnItsGene: StoryObj = { render: () => exampleStory(view) };
