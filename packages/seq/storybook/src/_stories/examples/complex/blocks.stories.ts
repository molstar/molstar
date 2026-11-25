/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';
import { bookmarkInterface, setCopies } from './compose.js';
import { complexEditStory } from './edit-root.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Complex/Blocks' };
export default meta;

/**
 * PDB 1FIN: CDK2 and cyclin A as blocks in a row, the polymers first and the two ATPs after them as one-position
 * blocks; the two copies of the complex merged, each block labeled ×2. Each feature is drawn on its own chain's block
 * only; the interface residues of both copies under them. Zoom into a block: the others leave the view.
 */
export const Row: StoryObj = { render: () => exampleStory(view('blocks')) };

/** The copies split: a block a chain and a ligand. */
export const Split: StoryObj = { render: () => exampleStory(view('blocks-split')) };

/** A block a row of the tracks, aligned at their first residue, one zoom for all. */
export const Stacked: StoryObj = { render: () => exampleStory(view('blocks-stacked')) };

/** Stacked with the copies split: a row a chain, the two ATP chips on one row of their own. */
export const StackedSplit: StoryObj = { render: () => exampleStory(view('blocks-stacked-split')) };

/** The blocks on a ring: a block an arc with a gap between, their labels outside with the features'. */
export const OnARing: StoryObj = { render: () => exampleStory(view('blocks-ring')) };

/**
 * The editing view: bookmark the interface (both sides' interface residues) as a ref, merge and split the copies, each
 * one patch; Undo and Redo walk them. Hover the merged block: every member lights.
 */
export const Editing: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('blocks'),
      edits: [
        { label: 'Bookmark the interface', run: bookmarkInterface },
        { label: 'Split the copies', run: (s) => setCopies(s, 'complex', 'each') },
        { label: 'Merge the copies', run: (s) => setCopies(s, 'complex', 'merge') },
      ],
    }),
};
