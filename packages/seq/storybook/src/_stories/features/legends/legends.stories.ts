/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { FeatureViews } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { LegendPlacementsRoot } from './legends-root.js';

// The legend part placed four ways by its host (the legends design): a track header's line (compact, the one track's), a
// side panel (the rest, stacked), a band under the canvas (compact, in a row) and the export's band below its figure.
// The same entries everywhere: hover or click a swatch to highlight or select its items.

const view = FeatureViews.find((v) => v.id === 'themes-editing')!;

const meta: Meta = { title: 'Features/Legends' };
export default meta;

/** The legend part placed four ways: a track header's line (compact), a side panel, a band under the canvas and the export's. */
export const LegendPlacements: StoryObj = {
  render: () => {
    const el = document.createElement(LegendPlacementsRoot.tag) as LegendPlacementsRoot;
    el.model = view;
    return el;
  },
};
