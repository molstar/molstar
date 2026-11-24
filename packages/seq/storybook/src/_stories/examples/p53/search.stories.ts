/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = P53Views.find((v) => v.id === 'protein-on-gene')!;

const meta: Meta = { title: 'Examples/p53/Search' };
export default meta;

/**
 * The protein on its gene with the stories' search box (the application's: the library ships the search's model and
 * commands): PXXP, p53's proline-rich motifs (IUPAC: X any residue), is searched at once, every hit tinted on the
 * protein and lit on the gene through the codon map (projected, not translated). Enter or Next steps to the next hit
 * (Shift+Enter back; Mod+G and F3 with the view focused), the count says where you are; Save keeps the hits as a set
 * (one undo step), Esc or Clear ends it. Try a DNA site on the gene (GGATCC, both strands by default) or a regular
 * expression (`P.{2}P`).
 */
export const Search: StoryObj = { render: () => exampleStory(view, { search: { initial: 'PXXP' } }) };
