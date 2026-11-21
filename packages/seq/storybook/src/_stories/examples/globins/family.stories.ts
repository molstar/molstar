/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GlobinsViews } from '@molstar/seqx-examples';
import { viewerItem } from '@molstar/seqx-viewer/config';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => GlobinsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Globins/Rows' };
export default meta;

/**
 * The Pfam globin seed's 60 rows by its tree (neighbor joining on identity), the conservation logo above them, each row's
 * species (UniProt's) and its identity to the reference row (wild yak's α) beside it. Hover a cell: its column lights in every row and the
 * row's label. Click a tree node to fold its clade into one summary row; the wheel over the rows scrolls them, `]` and
 * `[` zoom them.
 */
export const Family: StoryObj = { render: () => exampleStory(view('family')) };

/**
 * The family with its positions rail under the last track (`navigator.railPlace: bottom`), the other place beside the
 * default one (under the ruler): the rail is the precise horizontal control beside the vertical overview's picture;
 * drag either, the other follows. `navigator.positionsRail: unlessNavigator` takes it away again.
 */
export const FamilyRailAtBottom: StoryObj = {
  render: () => exampleStory(view('family'), { config: [[viewerItem('navigator.railPlace'), 'bottom']] as never }),
};

/** Two clades folded from the spec (`collapsed`): a summary row each, its consensus opaque by its share, a tree wedge. */
export const Folded: StoryObj = { render: () => exampleStory(view('folded')) };

/** The seed as a text figure: every row on every line, groups of ten, the ruler above. */
export const Figure: StoryObj = { render: () => exampleStory(view('figure')) };

/**
 * 10,000 held rows sampled from the seed's per-column profile, by a balanced tree: zoomed out the rows are image tiles
 * of the shared reduction.
 */
export const Scale: StoryObj = { render: () => exampleStory(view('scale')) };

/**
 * The 10,000 rows with a vertical overview of every row at the right (the shared reduction's image) and a miniature of
 * the columns and rows above: drag either window to move the rows (the miniature's both axes), click outside to centre
 * there, PageDown and PageUp over a part move between its windows. The SVG includes both, titled.
 */
export const ScaleWithNavigators: StoryObj = {
  render: () => exampleStory(view('scale-navigators'), { navigators: 'include' }),
};

/** 10,000,000 lazy rows: only the rows in view are made, a row not yet read labeled by its number. */
export const Lazy: StoryObj = { render: () => exampleStory(view('lazy')) };
