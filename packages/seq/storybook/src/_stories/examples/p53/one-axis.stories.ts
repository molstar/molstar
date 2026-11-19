/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'one-axis')!;

const meta: Meta = { title: 'Examples/p53/Protein and gene on one axis' };
export default meta;

/**
 * One gene section: the protein's letters a codon wide, its domains and the gene's letters, read along the protein on
 * its coding pieces. Each intron is a hairline with a notch (hover it: the hidden stretch); a click on a residue
 * selects it, its codon shown; Shift+R reverses the section, Shift+C shows the introns.
 */
export const ProteinAndGeneOnOneAxis: StoryObj = { render: () => exampleStory(view) };
