/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

// An alignment's rows as a feature, on the Globins family: the canvas and the SVG export of one spec side by side, in
// the light and the dark theme. Look at what the summary says.

const meta: Meta = { title: 'Features/Rows' };
export default meta;

/**
 * Data columns beside the rows of the Pfam globin seed: the species as text and the class as a chip on the left
 * (UniProt's, from a file joined by row key; an organism with no class rank has an empty cell), the identity to the
 * reference row, the residues as a bar and the UniProt entry's whole length on the right, after the last numbers. Hover a cut species: it reads whole. Click a title to order the rows by
 * it (an arrow marks it), again for the other direction; undo goes back. Narrow the window: the right columns yield
 * first, each down to a floor and then hidden, then the left; a chip above the rows counts those hidden.
 */
export const Columns: StoryObj = {
  render: () => exampleStory(FeatureViews.find((v) => v.id === 'columns')!),
};
