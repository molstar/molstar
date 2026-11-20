/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MoleculeViews } from '@molstar/seqx-examples';
import { LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';

// Molecules drawn (Design › Molecules): the canvas and the SVG of one spec side by side, in the light and the dark
// theme. `linksHover` is on in every story (a host's choice: off in the default viewer).

const meta: Meta = { title: 'Features/Molecules' };
export default meta;

/** The host's choice the stories show on: hovering a molecule lights its sites, a site its molecule. */
const LinksHoverOn = [{ provider: LinksHoverProvider, params: undefined }] as const;

const story = (id: string) => {
  const view = MoleculeViews.find((v) => v.id === id)!;
  return { render: () => exampleStory(view, { behaviors: LinksHoverOn }) } satisfies StoryObj;
};

/**
 * A section of a complex: the two Fc chains merged into one block “Fc A ×2”, then the glycan, ATP and the sodium ion as
 * chips (the water left out). MSE at 1 and SEP at 21 carry a dot at their cell's top right, CRO a bar over its three
 * letters. Hover the glycan's chip: N47 lights; hover N47: the glycan lights. The arcs: the glycan's, the disulfide, the
 * metal. Zoom in: the chips keep their width.
 */
export const SectionOfAComplex: StoryObj = story('molecules-complex');

/** The forms: the molecules whole, then the glycan's residues by their codes and ATP's atoms by their elements. */
export const WholeResiduesAtoms: StoryObj = story('molecules-forms');

/**
 * A modified residue four ways on a tRNA: a dot (2MG at 10, H2U | U at 17, PSU at 39 and 55), the code in place of the
 * letter where it fits (zoom in), the cell outlined, nothing. Zoom out: the dot becomes a tick at the strip's top. Hover
 * one: its component and name; H2U | U says both.
 */
export const ModifiedResidues: StoryObj = story('molecules-modified');

/** As text: each chain in lines, then the glycan, ATP and the ion as chips flowing like words on one line. */
export const AsText: StoryObj = story('molecules-text');

/** 200 copies of the chain as one block “Fc A ×200”, the molecules' chips after it: one block's cost. */
export const TwoHundredCopies: StoryObj = story('molecules-copies');
