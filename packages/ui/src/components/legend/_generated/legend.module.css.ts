/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from legend.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'bar': 'msui-legend-bar',
  'caption': 'msui-legend-caption',
  'dot': 'msui-legend-dot',
  'ends': 'msui-legend-ends',
  'item': 'msui-legend-item',
  'label': 'msui-legend-label',
  'legend': 'msui-legend-legend',
  'line': 'msui-legend-line',
  'mark': 'msui-legend-mark',
  'none': 'msui-legend-none',
  'ramp': 'msui-legend-ramp',
  'slot': 'msui-legend-slot',
  'swatch': 'msui-legend-swatch',
  'title': 'msui-legend-title',
} as const;

defineStyles('components/legend/legend.module.css', `@layer msui.components {
.msui-legend-legend {
  gap: var(--msui-space-1);
  min-width: 0;
  color: var(--msui-text);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
  flex-direction: column;
  display: flex;
}

.msui-legend-title {
  color: var(--msui-text-muted);
  font-weight: var(--msui-weight-strong);
}

.msui-legend-item {
  grid-template-columns: var(--msui-icon-slot) minmax(0, 1fr);
  column-gap: var(--msui-icon-gap);
  align-items: center;
  min-width: 0;
  display: grid;
}

.msui-legend-slot {
  width: var(--msui-icon-slot);
  height: var(--msui-icon);
  place-items: center;
  display: grid;
}

.msui-legend-label {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-legend-mark {
  background: var(--msui-legend-color, transparent);
}

.msui-legend-swatch {
  border-radius: var(--msui-radius-sm);
  width: 12px;
  height: 12px;
}

.msui-legend-dot {
  border-radius: 50%;
  width: 8px;
  height: 8px;
}

.msui-legend-line {
  width: var(--msui-icon);
  border-radius: var(--msui-radius-sm);
  height: 4px;
}

.msui-legend-none {
  background: color-mix(in oklab, currentColor, transparent 85%);
}

.msui-legend-ramp {
  gap: var(--msui-space-1);
  flex-direction: column;
  display: flex;
}

.msui-legend-caption {
  color: var(--msui-text-muted);
}

.msui-legend-bar {
  border-radius: var(--msui-radius-sm);
  height: 8px;
  display: block;
}

.msui-legend-ends {
  justify-content: space-between;
  gap: var(--msui-space-3);
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  font-variant-numeric: tabular-nums;
  display: flex;
}
}
`);
