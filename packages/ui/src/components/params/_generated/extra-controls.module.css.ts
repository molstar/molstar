/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from extra-controls.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'chevron': 'msui-params-extra-controls-chevron',
  'curve': 'msui-params-extra-controls-curve',
  'field': 'msui-params-extra-controls-field',
  'fieldContent': 'msui-params-extra-controls-fieldContent',
  'graph': 'msui-params-extra-controls-graph',
  'graphArea': 'msui-params-extra-controls-graphArea',
  'graphAt': 'msui-params-extra-controls-graphAt',
  'graphEditor': 'msui-params-extra-controls-graphEditor',
  'graphFoot': 'msui-params-extra-controls-graphFoot',
  'graphHint': 'msui-params-extra-controls-graphHint',
  'grid': 'msui-params-extra-controls-grid',
  'interval': 'msui-params-extra-controls-interval',
  'intervalRange': 'msui-params-extra-controls-intervalRange',
  'intervalTrack': 'msui-params-extra-controls-intervalTrack',
  'intervalValue': 'msui-params-extra-controls-intervalValue',
  'matrix': 'msui-params-extra-controls-matrix',
  'optionContent': 'msui-params-extra-controls-optionContent',
  'pair': 'msui-params-extra-controls-pair',
  'panel': 'msui-params-extra-controls-panel',
  'point': 'msui-params-extra-controls-point',
  'spark': 'msui-params-extra-controls-spark',
  'strip': 'msui-params-extra-controls-strip',
  'value': 'msui-params-extra-controls-value',
  'vec': 'msui-params-extra-controls-vec',
} as const;

defineStyles('components/params/extra-controls.module.css', `@layer msui.components {
.msui-params-extra-controls-field {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: center;
  gap: var(--msui-icon-gap);
  width: 100%;
  min-width: 0;
  height: var(--msui-control-height);
  padding: 0 var(--msui-space-1) 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  color: var(--msui-text);
  text-align: start;
  display: flex;

  &:hover, &[aria-expanded="true"] {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-params-extra-controls-fieldContent {
  align-items: center;
  gap: var(--msui-space-2);
  flex: 1;
  min-width: 0;
  display: flex;
}

.msui-params-extra-controls-value {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  flex: 1;
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-params-extra-controls-chevron {
  color: var(--msui-text-muted);
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  flex: none;

  [aria-expanded="true"] > & {
    transform: rotate(180deg);
  }
}

.msui-params-extra-controls-panel {
  gap: var(--msui-space-2);
  padding: var(--msui-space-1) var(--msui-inner) 0;
  display: grid;
}

.msui-params-extra-controls-optionContent {
  align-items: center;
  gap: var(--msui-space-2);
  min-width: 0;
  display: flex;
}

.msui-params-extra-controls-vec {
  gap: var(--msui-space-1);
  grid-template-columns: repeat(3, minmax(0, 1fr));
  width: 100%;
  display: grid;
}

.msui-params-extra-controls-matrix {
  gap: var(--msui-space-1);
  grid-template-columns: repeat(4, minmax(0, 1fr));
  display: grid;
}

.msui-params-extra-controls-interval {
  align-items: center;
  gap: var(--msui-space-2);
  grid-template-columns: minmax(0, 1fr) 72px;
  width: 100%;
  display: grid;
}

.msui-params-extra-controls-intervalTrack {
  height: var(--msui-control-height);
  margin-inline-start: var(--msui-control-inset);
  position: relative;

  &:before {
    content: "";
    border-radius: var(--msui-radius-sm);
    background: linear-gradient(90deg,
      var(--msui-control-active) var(--msui-params-from),
      color-mix(in oklab, var(--msui-text-muted), transparent 45%) var(--msui-params-from) var(--msui-params-to),
      var(--msui-control-active) var(--msui-params-to));
    height: 4px;
    position: absolute;
    top: 50%;
    left: 0;
    right: 0;
    translate: 0 -50%;
  }
}

.msui-params-extra-controls-intervalRange {
  -webkit-appearance: none;
  appearance: none;
  pointer-events: none;
  background: none;
  width: 100%;
  margin: 0;
  position: absolute;
  inset: 0;

  &::-webkit-slider-runnable-track {
    background: none;
    height: 4px;
  }

  &::-webkit-slider-thumb {
    background: var(--msui-text-muted);
    -webkit-appearance: none;
    appearance: none;
    pointer-events: auto;
    cursor: pointer;
    border-radius: 50%;
    width: 12px;
    height: 12px;
    margin-top: -4px;
  }

  &::-moz-range-track {
    background: none;
  }

  &::-moz-range-thumb {
    background: var(--msui-text-muted);
    pointer-events: auto;
    cursor: pointer;
    border: 0;
    border-radius: 50%;
    width: 12px;
    height: 12px;
  }

  &:focus-visible {
    outline: none;
  }

  &:focus-visible::-webkit-slider-thumb {
    outline: 2px solid var(--msui-focus);
    outline-offset: 2px;
  }
}

.msui-params-extra-controls-intervalValue {
  height: var(--msui-control-height);
  padding: 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  color: var(--msui-text);
  font-variant-numeric: tabular-nums;
  line-height: var(--msui-control-height);
  text-align: end;
  white-space: nowrap;
}

.msui-params-extra-controls-strip {
  border-radius: var(--msui-radius-sm);
  width: 48px;
  height: 8px;
  box-shadow: inset 0 0 0 1px var(--msui-line);
  flex: none;
}

.msui-params-extra-controls-spark {
  fill: none;
  stroke: currentColor;
  stroke-width: 1.2px;
  vector-effect: non-scaling-stroke;
  flex: none;
  width: 48px;
  height: 12px;
}

.msui-params-extra-controls-graphEditor {
  gap: var(--msui-space-1);
  display: grid;
}

.msui-params-extra-controls-graphArea {
  position: relative;
}

.msui-params-extra-controls-graph {
  border-radius: var(--msui-radius-md);
  background: var(--msui-surface-sunken);
  cursor: crosshair;
  touch-action: none;
  width: 100%;
  height: 96px;
  display: block;
}

.msui-params-extra-controls-grid {
  stroke: var(--msui-line);
  stroke-width: 1px;
  vector-effect: non-scaling-stroke;
}

.msui-params-extra-controls-curve {
  fill: none;
  stroke: var(--msui-text-muted);
  stroke-width: 1.5px;
  vector-effect: non-scaling-stroke;
}

.msui-params-extra-controls-point {
  fill: var(--msui-surface);
  stroke: var(--msui-text);
  stroke-width: 1.5px;
  vector-effect: non-scaling-stroke;
  cursor: grab;

  &:hover {
    fill: var(--msui-accent);
  }
}

.msui-params-extra-controls-graphAt {
  top: var(--msui-space-1);
  right: var(--msui-space-2);
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  font-variant-numeric: tabular-nums;
  pointer-events: none;
  position: absolute;
}

.msui-params-extra-controls-graphFoot {
  align-items: start;
  gap: var(--msui-space-2);
  min-width: 0;
  display: flex;
}

.msui-params-extra-controls-graphHint {
  min-width: 0;
  padding-block: var(--msui-space-1) 0;
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  flex: 1;
  margin: 0;
}

.msui-params-extra-controls-pair {
  gap: var(--msui-space-1);
  grid-template-columns: repeat(2, minmax(0, 1fr));
  width: 100%;
  display: grid;
}
}
`);
