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
 * The construct (a stand-in of pUC19's length) as rings, the plus strand's features outermost, then the minus
 * strand's, the sites and a ruler; labels outside with leaders, the name and length in the hole. The linear detail under
 * it follows the ring's window (the sector on the rim): wheel or pan the detail and the window moves round the ring,
 * across the origin too, where the detail shows two pieces. bla and a primer site cross the origin, each one arc.
 */
export const Ring: StoryObj = { render: () => exampleStory(view('ring')) };

/** The rotated linear view: tracks from inside bla over the origin to past lacZα, the seam no break. */
export const RotatedLinear: StoryObj = { render: () => exampleStory(view('rotated-linear')) };

/**
 * The ring and its linear detail with a strip of the whole construct over the detail: pan the detail over the origin
 * and the strip's window wraps in two pieces, a grip at each true end; drag either piece and the detail follows.
 */
export const DetailWithAStrip: StoryObj = {
  render: () => exampleStory(view('detail-strip'), { navigators: 'include' }),
};
