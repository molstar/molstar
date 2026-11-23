/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { RegionReadsViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { exampleStory } from '../../examples/example-root.js';
import { openEmpty } from './openings.js';

const view = (id: string) => RegionReadsViews.find((v) => v.id === id)!;

const meta: Meta = { title: 'Features/Region reads/Zoom in' };
export default meta;

/**
 * Two bands the view opens on at 400 to 700 kb. The genotypes (served at their items’ level alone) are past what a view
 * asks the store for even at the coarsest level: the band is hatched and says “Zoom in to read: ⌘/Ctrl+scroll, or +” (an info issue
 * in the status line), nothing is read and no wait resolves it; zoom in to a stretch of about 20 kb and it is read. The
 * links have no end in the stretch in view: nothing is drawn, no arc and no placeholder (a bin with nothing in it is not
 * data waiting). Pan left to the first 200 kb and the arcs are there (EC-169, EC-208).
 */
export const ZoomInToRead: StoryObj = {
  render: () => exampleStory(view('zoom-in'), { stacked: true, opening: openEmpty }),
};
