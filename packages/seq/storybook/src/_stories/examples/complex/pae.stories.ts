/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ComplexViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => ComplexViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Complex/PAE and contacts' };
export default meta;

/**
 * A section a chain (CDK2, cyclin A; stand-in letters at 1FIN's lengths) and the predicted aligned error between them
 * as a grid: its x follows CDK2's section, its y cyclin A's, so zooming either chain zooms that axis of the plot. The
 * interface's contacts are dots on it. Hover a cell: a crosshair, and both residues lit in their sections.
 */
export const Pae: StoryObj = { render: () => exampleStory(view('pae')) };

/** The interface's contacts as arcs from CDK2's block to cyclin A's: hover one, both ends light. */
export const Contacts: StoryObj = { render: () => exampleStory(view('contacts')) };
