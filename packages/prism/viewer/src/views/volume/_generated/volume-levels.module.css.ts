/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from volume-levels.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'axis': 'msp-viewer-volume-levels-axis',
  'bars': 'msp-viewer-volume-levels-bars',
  'dim': 'msp-viewer-volume-levels-dim',
  'editor': 'msp-viewer-volume-levels-editor',
  'head': 'msp-viewer-volume-levels-head',
  'inside': 'msp-viewer-volume-levels-inside',
  'level': 'msp-viewer-volume-levels-level',
  'levels': 'msp-viewer-volume-levels-levels',
  'marker': 'msp-viewer-volume-levels-marker',
  'mean': 'msp-viewer-volume-levels-mean',
  'name': 'msp-viewer-volume-levels-name',
  'off': 'msp-viewer-volume-levels-off',
  'other': 'msp-viewer-volume-levels-other',
  'plot': 'msp-viewer-volume-levels-plot',
  'row': 'msp-viewer-volume-levels-row',
  'suffix': 'msp-viewer-volume-levels-suffix',
  'swatch': 'msp-viewer-volume-levels-swatch',
  'themed': 'msp-viewer-volume-levels-themed',
  'title': 'msp-viewer-volume-levels-title',
  'track': 'msp-viewer-volume-levels-track',
  'unit': 'msp-viewer-volume-levels-unit',
  'units': 'msp-viewer-volume-levels-units',
  'value': 'msp-viewer-volume-levels-value',
} as const;

defineStyles('viewer/views/volume/volume-levels.module.css', `@layer msp.components {
.msp-viewer-volume-levels-levels {
  gap: var(--msp-space-1);
  font-size: var(--msp-font-size-sm);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-volume-levels-head {
  justify-content: space-between;
  align-items: center;
  gap: var(--msp-space-2);
  display: flex;
}

.msp-viewer-volume-levels-title {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-weight: 500;
}

.msp-viewer-volume-levels-units {
  gap: var(--msp-space-1);
  display: flex;
}

.msp-viewer-volume-levels-unit {
  min-width: 3ch;
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  line-height: 18px;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-checked="true"] {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-viewer-volume-levels-plot {
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-well);
  cursor: ew-resize;
  touch-action: none;
  contain: strict;
  height: 48px;
  position: relative;
  overflow: hidden;
}

.msp-viewer-volume-levels-bars {
  width: 100%;
  height: 100%;
  position: absolute;
  inset: 0;

  & path {
    fill: color-mix(in oklab, var(--msp-color-text-muted), transparent 60%);
  }
}

.msp-viewer-volume-levels-mean {
  pointer-events: none;
  position: absolute;
  inset: 0;

  &:before {
    content: "";
    border-left: 1px dashed color-mix(in oklab, var(--msp-color-text-subtle), transparent 40%);
    position: absolute;
    top: 0;
    bottom: 0;
    left: 0;
  }
}

.msp-viewer-volume-levels-inside {
  background: color-mix(in oklab, var(--msp-level, var(--msp-color-text-muted)), transparent 82%);
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msp-viewer-volume-levels-track {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msp-viewer-volume-levels-marker {
  pointer-events: auto;
  cursor: ew-resize;
  outline: none;
  width: 10px;
  position: absolute;
  top: 0;
  bottom: 0;
  left: -5px;

  &:before {
    content: "";
    background: var(--msp-level, var(--msp-color-text));
    width: 2px;
    position: absolute;
    top: 0;
    bottom: 0;
    left: 4px;
  }

  &:after {
    content: "";
    background: var(--msp-level, var(--msp-color-text));
    width: 8px;
    height: 8px;
    box-shadow: 0 0 0 1.5px var(--msp-color-well);
    transition: transform var(--msp-motion-fast) var(--msp-ease);
    border-radius: 50%;
    position: absolute;
    top: 2px;
    left: 1px;
  }

  &:hover:after, &:focus-visible:after {
    transform: scale(1.3);
  }

  &:focus-visible:before {
    box-shadow: 0 0 0 1px var(--msp-color-focus-ring);
  }
}

.msp-viewer-volume-levels-off {
  opacity: .35;
}

.msp-viewer-volume-levels-axis {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  justify-content: space-between;
  display: flex;
}

.msp-viewer-volume-levels-level {
  display: grid;
}

.msp-viewer-volume-levels-row {
  align-items: center;
  gap: var(--msp-space-2);
  height: var(--msp-row-height);
  display: flex;
}

.msp-viewer-volume-levels-swatch {
  border-radius: var(--msp-radius-sm);
  background: var(--msp-level, var(--msp-color-text-muted));
  width: 14px;
  height: 14px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 85%);
  flex: none;

  &:hover {
    transform: scale(1.1);
  }
}

.msp-viewer-volume-levels-themed {
  background: conic-gradient(var(--msp-accent-red) 0 25%,
    var(--msp-accent-orange) 0 50%,
    var(--msp-accent-green) 0 75%,
    var(--msp-accent-blue) 0);
}

.msp-viewer-volume-levels-name {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-volume-levels-dim {
  color: var(--msp-color-text-subtle);
}

.msp-viewer-volume-levels-value {
  width: 7ch;
  height: 22px;
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  color: inherit;
  font: inherit;
  font-variant-numeric: tabular-nums;
  text-align: right;
  border: 0;
  flex: none;

  &:focus-visible {
    outline: 1px solid var(--msp-color-focus-ring);
  }
}

.msp-viewer-volume-levels-suffix {
  width: 1ch;
  color: var(--msp-color-text-subtle);
  flex: none;
}

.msp-viewer-volume-levels-other {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  text-align: right;
  flex: none;
}

.msp-viewer-volume-levels-editor {
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) 0 var(--msp-space-2);
  display: grid;
}
}
`);
