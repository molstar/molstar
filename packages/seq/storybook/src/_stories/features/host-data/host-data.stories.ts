/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { hostDataStory } from './host-data-root.js';

// Host data (Design/Seq data ports): a spec that names PDB 2OCJ by its file, read from a host's memory instead.

const meta: Meta = { title: 'Features/Host data' };
export default meta;

/**
 * The fake structure viewer holds 2OCJ and Seq* is given its adapter: the nodes (the sequence, the zinc, the compounds,
 * the complex with its eight units) and the B-factors of chains A and D come from the fake's chain table; "What Seq*
 * holds" says that no asset was read. "Drop chain" (the last chain, D, or the chain of the cell you selected) edits the
 * fake's table: the source reopens a frame later and the listing follows, the other chains' nodes kept; chain D's
 * track empties and the listing under the view says why (its unit is gone from the source: `seq.resolve.system.unnamed`, and the
 * section has no length to draw: `seq.render.axis.pending` at `/views/main/chain-D`). The button above reads the file instead, with the default adapter (the fake starts again from
 * its table): the same spec, the same nodes, one asset read.
 */
export const ServedByTheFake: StoryObj = {
  render: () =>
    hostDataStory({
      hint: 'Seq* shows the B-factors of chains A and D of PDB 2OCJ, read from the fake structure viewer’s memory: no file is read. The legend under the view says the bars’ scale (Å²) and why the tracks end short of the sequence (the data). Click "Drop chain" in the fake: the source reopens a frame later; the listing under it follows, and chain D’s track empties with an issue that says its unit is gone. "Read the file instead" composes the same spec with the default adapter and starts the fake again.',
    }),
};

/**
 * The host's own complex as a section: the spec names the complex the source holds (its dynamic node), and the section
 * shows a block a unit, chains A to D and their zinc ions, once the source opens (pending before, never an error).
 * "Drop chain" in the fake removes a chain and its zinc; the section reflows a frame later, the other blocks taking the
 * room; "Restore the chains" brings them back the same way. The listing under the view prints the section's blocks.
 */
export const TheHostsComplex: StoryObj = {
  render: () =>
    hostDataStory({
      spec: 'complex',
      hint: 'Seq* shows the complex the fake structure viewer holds (PDB 2OCJ) as one section: a block a unit, the four chains and their zinc ions. Click "Drop chain" in the fake: a frame later the chain’s blocks leave and the others take the room. "Restore the chains" brings them back. No file is read.',
    }),
};
