/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { GlobinsViews } from '@molstar/seqx-examples';
import { LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { complexEditStory } from '../complex/edit-root.js';
import { exampleStory } from '../example-root.js';
import { toggleCopies } from './editing.js';

const view = (id: string) => GlobinsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Globins/Hemoglobin' };
export default meta;

/** The host's choice the stories show on (off in the default viewer): a heme lights its contacts, a contact its heme. */
const LinksHoverOn = [{ provider: LinksHoverProvider, params: undefined }] as const;

/**
 * PDB 4HHB as a complex's units in blocks: the α chains A and C and the β chains B and D, each pair merged into one
 * block labeled ×2, the four hemes (the heme molecule, HEM with its 43 heavy atoms) one chip ×4 after them. The helices A to H from the entry; each heme's contacting
 * residues (within 4 Å, derived per pair of units) as arcs and as the heme pocket, bookmarked (the bookmarks list runs
 * each). Hover α: its copy is the same block; hover a contact's arc: the heme and its residue light; hover the heme's
 * chip: the residues it contacts light (`links.hover`, on here), and His87 lights the heme. The toggle splits
 * the copies into a block a chain, each heme's arcs to the chain it sits in, and merges them again (one patch each way,
 * undone and redone from the stack beside it).
 */
export const Hemoglobin: StoryObj = {
  render: () =>
    complexEditStory({
      view: view('hemoglobin'),
      edits: [{ label: 'Split/merge the copies', run: (s) => toggleCopies(s) }],
      behaviors: LinksHoverOn,
    }),
};

/** The copies split (`copies: each`): a block a chain, a chip a heme. */
export const CopiesSplit: StoryObj = {
  render: () => exampleStory(view('hemoglobin-split'), { behaviors: LinksHoverOn }),
};
