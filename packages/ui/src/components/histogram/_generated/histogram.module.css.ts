/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from histogram.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'axis': 'msui-histogram-axis',
  'bars': 'msui-histogram-bars',
  'control': 'msui-histogram-control',
  'head': 'msui-histogram-head',
  'layer': 'msui-histogram-layer',
  'mean': 'msui-histogram-mean',
  'plot': 'msui-histogram-plot',
  'title': 'msui-histogram-title',
  'tools': 'msui-histogram-tools',
  'unit': 'msui-histogram-unit',
  'units': 'msui-histogram-units',
} as const;

defineStyles('components/histogram/histogram.module.css', `@layer msui.components {
.msui-histogram-control {
  padding-inline: var(--msui-inner);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msui-histogram-head {
  justify-content: space-between;
  align-items: center;
  gap: var(--msui-space-2);
  height: var(--msui-row);
  display: flex;
}

.msui-histogram-title {
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  font-weight: var(--msui-weight-strong);
  white-space: nowrap;
}

.msui-histogram-tools {
  align-items: center;
  gap: var(--msui-space-1);
  display: flex;
}

.msui-histogram-units {
  gap: var(--msui-space-1);
  display: flex;
}

.msui-histogram-unit {
  min-width: var(--msui-control-height);
  height: var(--msui-control-height);
  padding: 0 var(--msui-space-1);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  white-space: nowrap;
  transition: background var(--msui-motion-fast) var(--msui-ease),
    color var(--msui-motion-fast) var(--msui-ease);
  justify-content: center;
  align-items: center;
  display: inline-flex;

  &:hover {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &[aria-checked="true"] {
    background: var(--msui-control-active);
    color: var(--msui-text);
    box-shadow: inset 0 0 0 1px var(--msui-selected-edge);
  }

  &:focus-visible {
    outline: 1px solid var(--msui-focus);
  }
}

.msui-histogram-plot {
  height: var(--msui-plot-height, 48px);
  border-radius: var(--msui-radius-sm);
  background: var(--msui-surface-sunken);
  touch-action: none;
  contain: strict;
  position: relative;
  overflow: hidden;
}

.msui-histogram-layer {
  width: 100%;
  height: 100%;
  position: absolute;
  inset: 0;
  overflow: visible;
}

.msui-histogram-bars {
  fill: color-mix(in oklab, var(--msui-text-muted), transparent 70%);
}

.msui-histogram-mean {
  pointer-events: none;
  position: absolute;
  inset: 0;

  &:before {
    content: "";
    border-left: 1px dashed color-mix(in oklab, var(--msui-text-subtle), transparent 30%);
    position: absolute;
    top: 0;
    bottom: 0;
    left: 0;
  }
}

.msui-histogram-axis {
  justify-content: space-between;
  gap: var(--msui-space-2);
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  font-variant-numeric: tabular-nums;
  display: flex;
}
}
`);
