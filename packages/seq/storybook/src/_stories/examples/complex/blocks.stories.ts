/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Complex/Blocks' };
export default meta;

/**
 * CDK2 and cyclin A (stand-in letters at 1FIN's lengths) as blocks in a row, the polymers first and ATP after them as a
 * one-position block; the two copies of the complex merged, each block labeled ×2. Each feature is drawn on its own
 * chain's block only. Zoom into a block: the others leave the view.
 */
export const Row: StoryObj = { render: () => exampleStory(view('blocks')) };

/** The copies split: a block a chain and a ligand. */
export const Split: StoryObj = { render: () => exampleStory(view('blocks-split')) };

/** A block a row of the tracks, aligned at their first residue, one zoom for all. */
export const Stacked: StoryObj = { render: () => exampleStory(view('blocks-stacked')) };

/** The blocks on a ring: a block an arc with a gap between, their labels outside with the features'. */
export const OnARing: StoryObj = { render: () => exampleStory(view('blocks-ring')) };
