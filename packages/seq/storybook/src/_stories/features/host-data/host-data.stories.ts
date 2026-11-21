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
      hint: 'Seq* shows the B-factors of chains A and D of PDB 2OCJ, read from the fake structure viewer’s memory: no file is read. The legend under the view says the bars' scale (Å²) and why the tracks end short of the sequence (the data). Click "Drop chain" in the fake: the source reopens a frame later; the listing under it follows, and chain D’s track empties with an issue that says its unit is gone. "Read the file instead" composes the same spec with the default adapter and starts the fake again.',
    }),
};
