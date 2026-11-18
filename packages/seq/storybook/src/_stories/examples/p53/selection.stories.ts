/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { selectionStory } from '../../features/selection/selection-root.js';

// Selection on the protein on its gene: the Features selection story's root with this example's hint. A drag on the gene's
// letters selects nucleotides; the protein's section shows the same selection projected through the codon map, lighter,
// in the canvas and in the SVG beside it.

const meta: Meta = { title: 'Examples/p53/Selection' };
export default meta;

const Hint =
  'Drag across the gene’s letters (the lower section) to select nucleotides: the protein above shows the codons they fall in, lighter, and the SVG beside the canvas draws the same selection. Esc clears.';

/** A drag-select on the gene, projected onto the protein, in the canvas and the SVG. */
export const Selection: StoryObj = {
  render: () => selectionStory(P53Views.find((v) => v.id === 'protein-on-gene')!, Hint),
};
