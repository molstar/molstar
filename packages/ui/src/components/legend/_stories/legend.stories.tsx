/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Legend, type LegendItem, type LegendRamp } from '../legend.js';

const meta: Meta = { title: 'Components/Legend' };
export default meta;

const items: LegendItem[] = [
  { label: 'Helix', color: '#e1405f' },
  { label: 'Strand', color: '#f2b53a', glyph: 'dot' },
  { label: 'Coil', color: '#8a93a3', glyph: 'line' },
  { label: 'Not modeled', color: null },
];

const ramp: LegendRamp = {
  caption: 'Hydrophobicity',
  lo: 'Polar',
  hi: 'Hydrophobic',
  stops: [
    { t: 0, color: '#2c7bb6' },
    { t: 0.5, color: '#ffffbf' },
    { t: 1, color: '#d7191c' },
  ],
};

// The canvas carries the ink: a theme (`data-msui-theme`) and its base color, as a viewer sets them over its canvas.
const canvas = (background: string) =>
  `display: grid; gap: var(--msui-space-4); width: 240px; color: var(--msui-text); padding: var(--msui-space-4); background: ${background}; --msui-bg: ${background};`;

const On = (theme: 'light' | 'dark', background: string): StoryObj => ({
  render: () =>
    story(() => {
      const shown = state<readonly LegendItem[]>(items);
      return (
        <div data-msui-theme={theme} style={canvas(background)}>
          <Legend title="Secondary structure" items={shown} />
          <Legend ramps={[ramp]} />
        </div>
      );
    }),
});

/** Items (a swatch, a dot, a line, "no color") and a ramp over a light canvas. */
export const Light: StoryObj = On('light', '#f4f5f7');

/** The same over a dark canvas. */
export const Dark: StoryObj = On('dark', '#1b1d22');
