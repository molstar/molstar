/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GlobinsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => GlobinsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Globins/Rows' };
export default meta;

/**
 * The family's rows by its tree, the conservation logo above them. Hover a cell: its column lights in every row and the
 * row's label. Click a tree node to fold its clade into one summary row; the wheel over the rows scrolls them, `]` and
 * `[` zoom them.
 */
export const Family: StoryObj = { render: () => exampleStory(view('family')) };

/** Two clades folded from the spec (`collapsed`): a summary row each, its consensus opaque by its share, a tree wedge. */
export const Folded: StoryObj = { render: () => exampleStory(view('folded')) };

/** The family as a text figure: every row on every line, groups of ten, the ruler above. */
export const Figure: StoryObj = { render: () => exampleStory(view('figure')) };

/** 10,000 held rows by a balanced tree: zoomed out the rows are image tiles of the shared reduction. */
export const Scale: StoryObj = { render: () => exampleStory(view('scale')) };

/** 10,000,000 lazy rows: only the rows in view are made, a row not yet read labeled by its number. */
export const Lazy: StoryObj = { render: () => exampleStory(view('lazy')) };
