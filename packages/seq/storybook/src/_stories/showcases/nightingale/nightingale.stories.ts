/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { nightingaleStory } from './nightingale-root.js';

// Nightingale reading the resolved view (composition.md: other viewers read it): its web components, fed by an adapter
// from the resolved view and the store; the story's root and what it does are in `nightingale-root.ts`.

const meta: Meta = { title: 'Showcases/Nightingale' };
export default meta;

const tracks = () => P53Views.find((v) => v.id === 'tracks')!;

/**
 * p53's letters from the store's letters tile and its domains from the ranges tile, both found through the resolved
 * view's layers, fed to `nightingale-sequence` and `nightingale-track`; Seq*'s own export of the same positions beside it.
 */
export const SequenceAndFeatures: StoryObj = {
  render: () =>
    nightingaleStory({
      view: tracks(),
      kind: 'features',
      reads:
        'Read from the resolved view and the store: the letters layer’s tile gives the sequence, the ranges layer’s the domains (0-based, half-open here, 1-based and inclusive in Nightingale), colored by the layer’s theme form from Seq*’s tokens. Ctrl+wheel zooms Nightingale.',
    }),
};

/** The conservation values layer, read from the store's values tile and drawn by `nightingale-linegraph-track`. */
export const Values: StoryObj = {
  render: () =>
    nightingaleStory({
      view: tracks(),
      kind: 'values',
      reads:
        'Read from the store: the values layer’s tile gives one value a position (two positions have none and are left out of the line); the line takes the layer’s color. Seq*’s own export of the same data is beside it.',
    }),
};

/**
 * Seq*’s hover and selection drive Nightingale’s `highlight`: a location of the viewer is a range of the same
 * positions, so only the base (0 to 1) changes.
 */
export const TheLink: StoryObj = {
  render: () =>
    nightingaleStory({
      view: tracks(),
      kind: 'link',
      reads:
        'Hover or click a residue, a domain or a conservation bar in Seq*: the viewer’s hover and selection locations (an address on the p53 system) are read and set as the `highlight` of Nightingale’s sequence and track. One way only: Seq* to Nightingale.',
    }),
};
