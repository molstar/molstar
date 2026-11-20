/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import { LinksHoverProvider } from '@molstar/seqx-viewer/behaviors';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => P53Views.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/p53/Structure' };
export default meta;

/** The host's choice the stories show on (off in the default viewer): a zinc lights its ligands, a ligand its zinc. */
const LinksHoverOn = [{ provider: LinksHoverProvider, params: undefined }] as const;

/**
 * P04637 over PDB 2OCJ: the four copies of the DNA-binding domain merged into one block ×4, the zincs one chip ×4;
 * the ends (94–95, 290–312) unobserved and gray; the entity's helices and strands, each zinc's four ligands as arcs.
 * Hover R175 on P04637: it lights on the chains through SIFTS's map, and back. Hover the zinc's chip: C176, H179, C238
 * and C242 light (`links.hover`, on here); hover C176: the zinc lights.
 */
export const Structure: StoryObj = { render: () => exampleStory(view('structure'), { behaviors: LinksHoverOn }) };

/** The copies split: a block a chain, a chip a zinc; a hover carried through the alignment's map (`via`). */
export const CopiesSplit: StoryObj = {
  render: () => exampleStory(view('structure-split'), { behaviors: LinksHoverOn }),
};
