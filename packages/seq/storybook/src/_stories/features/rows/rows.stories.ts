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
 * Data columns beside the rows: the species as text and the class as a chip on the left (the class from a file joined
 * by row key, every fifth row without one: an empty cell), the identity to the ancestor, the residues as a bar and a
 * score on the right, after the last numbers. Hover a cut species: it reads whole. Click a title to order the rows by
 * it (an arrow marks it), again for the other direction; undo goes back. Narrow the window: the right columns yield
 * first, each down to a floor and then hidden, then the left; a chip above the rows counts those hidden.
 */
export const Columns: StoryObj = {
  render: () => exampleStory(FeatureViews.find((v) => v.id === 'columns')!),
};
