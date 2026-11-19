/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PlasmidViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../example-root.js';

const view = (id: string) => PlasmidViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Examples/Plasmid/Ring' };
export default meta;

/**
 * pUC19c (INSDC L09137.2, 2,686 bp) as rings: the plus strand's features outermost, then the minus strand's, the
 * restriction sites at their cuts and a ruler; labels outside with leaders, the name and length in the hole. Features
 * are arrows by strand (plus to the end, minus to the start; a ring draws them from a device pixel a position, so on a
 * high-resolution screen or large). The linear detail under it follows the ring's window (the sector on the rim):
 * wheel or pan the detail and the window moves round the ring, across the origin too, where the detail shows two
 * pieces. bla and a primer site cross the origin, each one arc.
 */
export const Ring: StoryObj = { render: () => exampleStory(view('ring')) };

/**
 * The ring and its linear detail with a strip of the whole construct over the detail: pan the detail over the origin
 * and the strip's window wraps in two pieces, a grip at each true end; drag either piece and the detail follows.
 */
export const DetailWithAStrip: StoryObj = {
  render: () => exampleStory(view('detail-strip'), { navigators: 'include' }),
};
