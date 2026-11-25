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
import { densityData, densityValues } from '../../histogram/_stories/density.js';
import { type ValueUnit, valueText } from '../../histogram/histogram.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { type IsoLevel, IsoLevelControl } from '../iso-level.js';

const meta: Meta = { title: 'Components/Iso level' };
export default meta;

const panel = 'width: 320px; margin: var(--msui-space-4);';

/** A difference map's values: symmetric around zero, with thin tails both ways (the + and − peaks). */
function differenceValues(n = 80_000): Float64Array {
  const base = densityValues(n, 0, 11);
  const tails = densityValues(n, 1, 13);
  return base.map((v, i) => (i % 40 === 0 ? (i % 80 === 0 ? 1 : -1) * (tails[i] - 0.1) * 0.6 : v));
}

/** An electron density map: one level at 1.5 σ on the long right tail. Drag on the plot or use the arrow keys. */
export const Density: StoryObj = {
  render: () =>
    story(() => {
      const { stats, histogram } = densityData();
      const levels = state<readonly IsoLevel[]>([
        { id: 'surface', label: 'Surface', value: stats.mean + 1.5 * stats.sigma, color: 0x3362b2 },
      ]);
      const unit = state<ValueUnit>('sigma');
      const commits = state(0);
      return (
        <div ref={(el: HTMLElement) => onCleanup(installTooltips(el))} class="msui-surface" style={panel}>
          <IsoLevelControl
            stats={stats}
            histogram={histogram}
            levels={levels}
            unit={unit}
            onChange={(id, value) => levels.set(levels.value.map((l) => (l.id === id ? { ...l, value } : l)))}
            onCommit={() => commits.set(commits.value + 1)}
          />
          <div class="msui-form">
            <span class="msui-label">Iso value</span>
            <span class="msui-value msui-num" style="text-align: start;">
              {derive([levels, unit], (ls, u) => valueText(u === 'sigma' ? 'absolute' : 'sigma', ls[0].value, stats))}
            </span>
            <span class="msui-label">Steps</span>
            <span class="msui-value">{derive(commits, String)}</span>
          </div>
        </div>
      );
    }),
};

/**
 * A difference map: a level each side of the mean in its color (what each encloses tinted: above the + level, below
 * the − one), the − level hidden; the unit shared with a second control.
 */
export const DifferenceMap: StoryObj = {
  render: () =>
    story(() => {
      const { stats, histogram } = densityData(differenceValues());
      const unit = state<ValueUnit>('sigma');
      const levels = state<readonly IsoLevel[]>([
        { id: 'plus', label: 'Fo-Fc (+)', value: stats.mean + 3 * stats.sigma, color: 0x2ca02c },
        { id: 'minus', label: 'Fo-Fc (−)', value: stats.mean - 3 * stats.sigma, color: 0xd62728, hidden: true },
      ]);
      const density = densityData();
      const second = state<readonly IsoLevel[]>([
        { id: 'surface', label: '2Fo-Fc', value: density.stats.mean + density.stats.sigma, color: 0x3362b2 },
      ]);
      const set = (s: typeof levels) => (id: string, value: number) =>
        s.set(s.value.map((l) => (l.id === id ? { ...l, value } : l)));
      return (
        <div ref={(el: HTMLElement) => onCleanup(installTooltips(el))} class="msui-surface" style={panel}>
          <IsoLevelControl
            title="Fo-Fc"
            stats={stats}
            histogram={histogram}
            levels={levels}
            unit={unit}
            onChange={set(levels)}
          />
          <IsoLevelControl
            title="2Fo-Fc"
            stats={density.stats}
            histogram={density.histogram}
            levels={second}
            unit={unit}
            onChange={set(second)}
          />
        </div>
      );
    }),
};

/** Before the values arrive (streamed density): an empty well; the marker (no color given: the text's) still moves. */
export const NoHistogramYet: StoryObj = {
  render: () =>
    story(() => {
      const { stats } = densityData();
      const levels = state<readonly IsoLevel[]>([{ id: 'surface', label: 'Surface', value: stats.mean + stats.sigma }]);
      return (
        <div class="msui-surface" style={panel}>
          <IsoLevelControl
            stats={stats}
            levels={levels}
            onChange={(id, value) => levels.set(levels.value.map((l) => (l.id === id ? { ...l, value } : l)))}
          />
        </div>
      );
    }),
};
