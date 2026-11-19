/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';
import { complexEditStory } from './edit-root.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Complex/PAE and contacts' };
export default meta;

/**
 * A section a chain of the prediction (CDK2 P24941, cyclin A2 P20248) and the synthetic predicted aligned error between
 * them as a grid: its x follows CDK2's section, its y cyclin A's, so zooming either chain zooms that axis of the plot.
 * 1FIN's contacts are dots on it. Hover a cell: a crosshair, and both residues lit in their sections.
 */
export const Pae: StoryObj = { render: () => exampleStory(view('pae')) };

/**
 * The PAE over the whole complex (synthetic, AlphaFold-shaped): the chains side by side as blocks, the plot following
 * them on both axes, the chain boundaries the blocks' gaps across it. Hover a cell: the viewer's card names both residues with
 * their chains and the value; drag a rectangle across a boundary: it is split at it, the card summarizes it while it
 * is dragged, and the panel (once also the hover's readout) says what is aligned and scored, with the mean PAE.
 */
export const OverTheComplex: StoryObj = {
  render: () => complexEditStory({ view: view('complex-grid'), edits: [], readout: true }),
};

/** 1FIN's contacts (≤ 4 Å, both copies) as arcs from CDK2's block to cyclin A's: hover one, both ends light. */
export const Contacts: StoryObj = { render: () => exampleStory(view('contacts')) };
