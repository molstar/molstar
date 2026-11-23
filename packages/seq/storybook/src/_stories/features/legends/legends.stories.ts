/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { createSeqSpec, theme } from '@molstar/seqx-core/lang';
import { type ExampleView, FeatureViews } from '@molstar/seqx-examples';
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

/** 40 values hashed to colors, listed as drawn: the legend lists 12 and a "+28" button, which opens the rest in place. */
const longLegend: ExampleView = {
  dataset: 'legends',
  id: 'long-legend',
  title: 'A long legend',
  summary: 'Forty chromosome-like values on one track: the legend lists twelve and a “+28” button opens the rest.',
  build() {
    const b = createSeqSpec({ title: 'A long legend' });
    const n = 40;
    const p = b.data.sequence('p', { alphabet: 'dna', letters: 'ACGT'.repeat(n * 5) });
    const d = b.data.ranges('d', {
      target: p,
      rows: {
        start: Array.from({ length: n }, (_, i) => i * 20),
        end: Array.from({ length: n }, (_, i) => i * 20 + 16),
        id: Array.from({ length: n }, (_, i) => `i${i}`),
        chromosome: Array.from({ length: n }, (_, i) => `chr${i + 1}`),
      },
    });
    b.views.view('main', (v) => {
      const s = v.section('p', { of: p });
      s.track('letters').layer('l', { type: 'letters' });
      s.track('chromosomes', { label: 'Chromosomes' }).layer('d', {
        type: 'ranges',
        data: d,
        color: theme.field('chromosome', { list: false }),
      });
    });
    return b.build();
  },
};

/** A legend over 40 values: its "+28" opens the full list in place (bounded, scrollable) and closes it again. */
export const LongLegend: StoryObj = {
  render: () => {
    const el = document.createElement(LegendPlacementsRoot.tag) as LegendPlacementsRoot;
    el.model = longLegend;
    return el;
  },
};
