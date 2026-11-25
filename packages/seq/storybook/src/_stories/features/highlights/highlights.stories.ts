/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { P53Views } from '@molstar/seqx-examples';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { highlightsStory } from './highlights-root.js';

// Highlights: overlapping channels stay apart (the story's root and what it does are in `highlights-root.ts`).

const meta: Meta = { title: 'Features/Highlights' };
export default meta;

/**
 * p53's tracks with three overlapping channels (a domain, a motif search result, a host's own) set and cleared from
 * buttons: where they overlap the band is split into stripes in the order they were set; each channel has a style, a
 * tint, a band along the foot or an outline; the legend lists the active ones, each removable.
 */
export const OverlappingChannels: StoryObj = {
  render: () => highlightsStory(P53Views.find((v) => v.id === 'tracks')!),
};
