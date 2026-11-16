/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from transfer-function.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'area': 'msui-transfer-function-area',
  'bar': 'msui-transfer-function-bar',
  'colors': 'msui-transfer-function-colors',
  'curve': 'msui-transfer-function-curve',
  'foot': 'msui-transfer-function-foot',
  'handle': 'msui-transfer-function-handle',
  'handleTrack': 'msui-transfer-function-handleTrack',
  'on': 'msui-transfer-function-on',
  'palette': 'msui-transfer-function-palette',
  'palettes': 'msui-transfer-function-palettes',
  'panel': 'msui-transfer-function-panel',
  'plot': 'msui-transfer-function-plot',
  'point': 'msui-transfer-function-point',
  'svg': 'msui-transfer-function-svg',
  'swatch': 'msui-transfer-function-swatch',
  'swatches': 'msui-transfer-function-swatches',
  'track': 'msui-transfer-function-track',
} as const;

defineStyles('components/transfer-function/transfer-function.module.css', `@layer msui.components {
.msui-transfer-function-plot {
  border-radius: var(--msui-radius-sm) var(--msui-radius-sm) 0 0;
  cursor: crosshair;
}

.msui-transfer-function-svg {
  width: 100%;
  height: 100%;
  position: absolute;
  inset: 0;
  overflow: visible;
}

.msui-transfer-function-area {
  opacity: .72;
}

.msui-transfer-function-curve {
  fill: none;
  stroke: var(--msui-text);
  stroke-width: 1.5px;
  stroke-linejoin: round;
  vector-effect: non-scaling-stroke;
}

.msui-transfer-function-track {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msui-transfer-function-point {
  top: calc(-1 * var(--msui-space-3) / 2);
  left: calc(-1 * var(--msui-space-3) / 2);
  width: var(--msui-space-3);
  height: var(--msui-space-3);
  pointer-events: auto;
  cursor: grab;
  outline: none;
  position: absolute;

  &:after {
    content: "";
    background: var(--msui-surface);
    box-shadow: inset 0 0 0 1.5px var(--msui-text);
    transition: transform var(--msui-motion-fast) var(--msui-ease);
    border-radius: 50%;
    position: absolute;
    inset: 2px;
  }

  &:hover:after, &:focus-visible:after {
    transform: scale(1.25);
  }

  &:focus-visible:after {
    box-shadow: inset 0 0 0 1.5px var(--msui-text),
      0 0 0 2px var(--msui-focus);
  }
}

.msui-transfer-function-on:after {
  background: var(--msui-accent);
  box-shadow: inset 0 0 0 1.5px var(--msui-accent);
}

.msui-transfer-function-colors {
  height: var(--msui-space-4);
  position: relative;
}

.msui-transfer-function-bar {
  height: var(--msui-space-2);
  border-radius: 0 0 var(--msui-radius-sm) var(--msui-radius-sm);
}

.msui-transfer-function-handleTrack {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msui-transfer-function-handle {
  top: var(--msui-space-1);
  left: calc(-1 * var(--msui-space-2) / 2);
  width: var(--msui-space-2);
  height: var(--msui-space-3);
  border-radius: var(--msui-radius-sm);
  background: var(--msui-stop, var(--msui-text-muted));
  box-shadow: 0 0 0 1px var(--msui-surface),
    0 0 0 2px color-mix(in oklab, var(--msui-text), transparent 60%);
  pointer-events: auto;
  cursor: pointer;
  border: 0;
  outline: none;
  padding: 0;
  position: absolute;

  &:hover {
    box-shadow: 0 0 0 1px var(--msui-surface),
      0 0 0 2px var(--msui-text);
  }

  &:focus-visible, &[aria-expanded="true"] {
    box-shadow: 0 0 0 1px var(--msui-surface),
      0 0 0 2px var(--msui-focus);
  }
}

.msui-transfer-function-palettes {
  display: grid;
}

.msui-transfer-function-palette {
  align-items: center;
  gap: var(--msui-space-2);
  height: var(--msui-row);
  padding-inline: var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text);
  text-align: start;
  white-space: nowrap;
  display: flex;

  &:hover, &:focus-visible {
    background: var(--msui-control-hover);
    outline: none;
  }
}

.msui-transfer-function-swatches {
  width: calc(16 * var(--msui-unit));
  height: var(--msui-space-3);
  border-radius: var(--msui-radius-sm);
  box-shadow: inset 0 0 0 1px var(--msui-line);
  flex: none;
  display: flex;
  overflow: hidden;
}

.msui-transfer-function-swatch {
  flex: 1;
}

.msui-transfer-function-panel {
  gap: var(--msui-space-2);
  padding: var(--msui-space-1) var(--msui-inner) 0;
  display: grid;
}

.msui-transfer-function-foot {
  justify-content: flex-end;
  display: flex;
}
}
`);
