/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { derive } from '@molstar/core/reactive';
import { onCleanup } from '../../../dom/index.js';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { densityData } from '../../histogram/_stories/density.js';
import type { ValueStats } from '../../histogram/histogram.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import {
  TransferFunction,
  type TransferPalette,
  type TransferPoint,
  type TransferPreset,
} from '../transfer-function.js';

const meta: Meta = { title: 'Components/Transfer function' };
export default meta;

const panel = 'width: 320px; margin: var(--msui-space-4);';

/** A peak of opacity from σ `a` to σ `b` (its top `alpha`), as curves made from a map's statistics are. */
function peak(stats: ValueStats, a: number, b: number, alpha: number) {
  const at = (s: number) => stats.mean + s * stats.sigma;
  const mid = (a + b) / 2;
  return [
    { x: at(a), alpha: 0 },
    { x: at(a + (mid - a) / 2), alpha },
    { x: at(b - (b - mid) / 2), alpha },
    { x: at(b), alpha: 0 },
  ];
}

const Palettes: readonly TransferPalette[] = [
  { label: 'Categorical', colors: [0x1f77b4, 0xff7f0e, 0x2ca02c, 0xd62728, 0x9467bd] },
  { label: 'Viridis', colors: [0x440154, 0x3b528b, 0x21918c, 0x5ec962, 0xfde725], sequential: true },
  { label: 'Red to blue', colors: [0xd62728, 0xf7f7f7, 0x1f77b4], sequential: true },
];

/**
 * A direct volume of an electron density map: the opacity (a log axis) over the histogram's long right tail, two
 * peaks in their colors. Drag a point, press on the plot to add one, double-click or Delete to remove one; a handle
 * under the bar opens its point's color; Presets and the palette recolor or reshape the curve.
 */
export const DirectVolume: StoryObj = {
  render: () =>
    story(() => {
      const { stats, histogram } = densityData();
      const presets: TransferPreset[] = [
        { label: 'Above 1 σ', points: [...peak(stats, 0.8, 2, 0.05), ...peak(stats, 2.2, 6, 0.2)] },
        { label: 'Narrow at 2 σ', points: peak(stats, 1.6, 2.4, 0.1) },
        {
          label: 'Ramp',
          points: [
            { x: stats.mean, alpha: 0 },
            { x: stats.max, alpha: 1 },
          ],
        },
      ];
      const points = state<readonly TransferPoint[]>([
        ...peak(stats, 0.8, 2, 0.05).map((p) => ({ ...p, color: 0x3362b2 })),
        ...peak(stats, 2.2, 6, 0.2).map((p) => ({ ...p, color: 0xff7f0e })),
      ]);
      const commits = state(0);
      return (
        <div ref={(el: HTMLElement) => onCleanup(installTooltips(el))} class="msui-surface" style={panel}>
          <TransferFunction
            stats={stats}
            histogram={histogram}
            points={points}
            logMin={1e-3}
            presets={presets}
            palettes={Palettes}
            onChange={(p) => points.set(p)}
            onCommit={() => commits.set(commits.value + 1)}
          />
          <div class="msui-form">
            <span class="msui-label">Points</span>
            <span class="msui-value">{derive(points, (p) => String(p.length))}</span>
            <span class="msui-label">Steps</span>
            <span class="msui-value">{derive(commits, String)}</span>
          </div>
        </div>
      );
    }),
};

/** A linear opacity axis (no `logMin`), a ramp in a sequential palette's colors; no presets or palettes. */
export const LinearRamp: StoryObj = {
  render: () =>
    story(() => {
      const { stats, histogram } = densityData();
      const points = state<readonly TransferPoint[]>([
        { x: stats.mean, alpha: 0, color: 0x440154 },
        { x: stats.mean + 2 * stats.sigma, alpha: 0.4, color: 0x21918c },
        { x: stats.max, alpha: 0.9, color: 0xfde725 },
      ]);
      return (
        <div ref={(el: HTMLElement) => onCleanup(installTooltips(el))} class="msui-surface" style={panel}>
          <TransferFunction
            title="Opacity"
            stats={stats}
            histogram={histogram}
            points={points}
            onChange={(p) => points.set(p)}
          />
        </div>
      );
    }),
};
