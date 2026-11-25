/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { selectionStory } from './selection-root.js';

// Selection: the default grammar on the two p53 views (the interaction grammar); the story's root and what it does are
// in `selection-root.ts`.

const meta: Meta = { title: 'Features/Selection' };
export default meta;

const story = (id: string, hint?: string) => selectionStory(P53Views.find((v) => v.id === id)!, hint);

const PairHint =
  'Drag across the gene’s letters (the lower section) to select nucleotides, Cmd/Ctrl+drag to add another stretch, Alt/Option+drag to cut one out: the protein above shows the codons they fall in, lighter, and the SVG beside the canvas draws the same selection. Pan or zoom on either section (its ruler, the wheel, a pinch): both move. Esc clears.';

/** p53's tracks: select residues, domains and ranges; the cursor; the SVG beside it draws the same selection. */
export const Tracks: StoryObj = { render: () => story('tracks') };

/**
 * The protein on its gene: a drag on the gene's letters selects nucleotides, and the protein shows the codons they fall
 * in, lighter, in the canvas and in the SVG beside it; what is hovered or selected on one section shows on the other.
 */
export const ProteinOnItsGene: StoryObj = { render: () => story('protein-on-gene', PairHint) };
