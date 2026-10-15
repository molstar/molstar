/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from volume-transfer.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'area': 'msp-viewer-volume-transfer-area',
  'axis': 'msp-viewer-volume-transfer-axis',
  'bars': 'msp-viewer-volume-transfer-bars',
  'chevron': 'msp-viewer-volume-transfer-chevron',
  'colorRow': 'msp-viewer-volume-transfer-colorRow',
  'curve': 'msp-viewer-volume-transfer-curve',
  'head': 'msp-viewer-volume-transfer-head',
  'iconTool': 'msp-viewer-volume-transfer-iconTool',
  'on': 'msp-viewer-volume-transfer-on',
  'palette': 'msp-viewer-volume-transfer-palette',
  'paletteLabel': 'msp-viewer-volume-transfer-paletteLabel',
  'paletteSwatch': 'msp-viewer-volume-transfer-paletteSwatch',
  'paletteSwatches': 'msp-viewer-volume-transfer-paletteSwatches',
  'palettes': 'msp-viewer-volume-transfer-palettes',
  'picked': 'msp-viewer-volume-transfer-picked',
  'pickedText': 'msp-viewer-volume-transfer-pickedText',
  'plot': 'msp-viewer-volume-transfer-plot',
  'point': 'msp-viewer-volume-transfer-point',
  'stop': 'msp-viewer-volume-transfer-stop',
  'stopTrack': 'msp-viewer-volume-transfer-stopTrack',
  'strip': 'msp-viewer-volume-transfer-strip',
  'svg': 'msp-viewer-volume-transfer-svg',
  'swatch': 'msp-viewer-volume-transfer-swatch',
  'themed': 'msp-viewer-volume-transfer-themed',
  'title': 'msp-viewer-volume-transfer-title',
  'tool': 'msp-viewer-volume-transfer-tool',
  'tools': 'msp-viewer-volume-transfer-tools',
  'track': 'msp-viewer-volume-transfer-track',
  'transfer': 'msp-viewer-volume-transfer-transfer',
  'uncolored': 'msp-viewer-volume-transfer-uncolored',
  'unit': 'msp-viewer-volume-transfer-unit',
  'units': 'msp-viewer-volume-transfer-units',
} as const;

defineStyles('viewer/views/volume/volume-transfer.module.css', `@layer msp.components {
.msp-viewer-volume-transfer-transfer {
  gap: var(--msp-space-1);
  font-size: var(--msp-font-size-sm);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-volume-transfer-head {
  justify-content: space-between;
  align-items: center;
  gap: var(--msp-space-2);
  display: flex;
}

.msp-viewer-volume-transfer-title {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-weight: 500;
}

.msp-viewer-volume-transfer-tools {
  align-items: center;
  gap: var(--msp-space-1);
  display: flex;
}

.msp-viewer-volume-transfer-tool {
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
  flex: none;
  line-height: 20px;

  &:hover:not(:disabled) {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[aria-expanded="true"] {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }

  &:disabled {
    opacity: .45;
  }
}

.msp-viewer-volume-transfer-units {
  gap: var(--msp-space-1);
  display: flex;
}

.msp-viewer-volume-transfer-unit {
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

.msp-viewer-volume-transfer-plot {
  border-radius: var(--msp-radius-sm) var(--msp-radius-sm) 0 0;
  background: var(--msp-color-well);
  cursor: crosshair;
  touch-action: none;
  contain: strict;
  height: 64px;
  position: relative;
  overflow: hidden;
}

.msp-viewer-volume-transfer-svg {
  width: 100%;
  height: 100%;
  position: absolute;
  inset: 0;
}

.msp-viewer-volume-transfer-bars {
  fill: color-mix(in oklab, var(--msp-color-text-muted), transparent 72%);
}

.msp-viewer-volume-transfer-area {
  opacity: .7;
}

.msp-viewer-volume-transfer-uncolored {
  fill: color-mix(in oklab, var(--msp-color-text-muted), transparent 50%);
}

.msp-viewer-volume-transfer-curve {
  fill: none;
  stroke: var(--msp-color-text);
  stroke-width: 1.5px;
  vector-effect: non-scaling-stroke;
}

.msp-viewer-volume-transfer-track {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msp-viewer-volume-transfer-point {
  pointer-events: auto;
  cursor: grab;
  outline: none;
  width: 12px;
  height: 12px;
  position: absolute;
  top: -6px;
  left: -6px;

  &:after {
    content: "";
    background: var(--msp-color-bg);
    box-shadow: inset 0 0 0 1.5px var(--msp-color-text);
    transition: transform var(--msp-motion-fast) var(--msp-ease);
    border-radius: 50%;
    position: absolute;
    inset: 2px;
  }

  &:hover:after, &:focus-visible:after {
    transform: scale(1.25);
  }
}

.msp-viewer-volume-transfer-on:after {
  background: var(--msp-color-accent);
  box-shadow: inset 0 0 0 1.5px var(--msp-color-accent);
}

.msp-viewer-volume-transfer-strip {
  border-radius: 0 0 var(--msp-radius-sm) var(--msp-radius-sm);
  cursor: copy;
  touch-action: none;
  height: 10px;
  position: relative;
}

.msp-viewer-volume-transfer-themed {
  background: repeating-linear-gradient(135deg,
    var(--msp-color-well) 0 4px,
    color-mix(in oklab, var(--msp-color-text-muted), transparent 70%) 4px 8px);
  cursor: default;
}

.msp-viewer-volume-transfer-stopTrack {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msp-viewer-volume-transfer-stop {
  background: var(--msp-stop);
  width: 9px;
  height: 12px;
  box-shadow: 0 0 0 1.5px var(--msp-color-bg),
    0 0 0 2.5px color-mix(in oklab, var(--msp-color-text), transparent 60%);
  pointer-events: auto;
  cursor: ew-resize;
  border-radius: 3px;
  outline: none;
  position: absolute;
  top: 2px;
  left: -4.5px;

  &:focus-visible {
    box-shadow: 0 0 0 1.5px var(--msp-color-bg),
      0 0 0 2.5px var(--msp-color-focus-ring);
  }

  &.msp-viewer-volume-transfer-on {
    box-shadow: 0 0 0 1.5px var(--msp-color-bg),
      0 0 0 3px var(--msp-color-accent);
  }

  &:after {
    content: none;
  }
}

.msp-viewer-volume-transfer-axis {
  margin-top: var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  justify-content: space-between;
  display: flex;
}

.msp-viewer-volume-transfer-iconTool {
  border-radius: var(--msp-radius-sm);
  width: 20px;
  height: 20px;
  color: var(--msp-color-text-muted);
  place-items: center;
  display: grid;

  &:hover:not(:disabled) {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &:disabled {
    opacity: .45;
  }
}

.msp-viewer-volume-transfer-picked {
  align-items: center;
  gap: var(--msp-space-1);
  min-height: 24px;
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-viewer-volume-transfer-pickedText {
  font-variant-numeric: tabular-nums;
  text-align: left;
  white-space: nowrap;
  flex: 1;
  min-width: 0;
  overflow: hidden;
}

.msp-viewer-volume-transfer-colorRow {
  align-items: center;
  gap: var(--msp-space-2);
  min-width: 0;
  height: 24px;
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  flex: 1;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-volume-transfer-swatch {
  border-radius: var(--msp-radius-sm);
  width: 12px;
  height: 12px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 85%);
  flex: none;
}

.msp-viewer-volume-transfer-chevron {
  color: var(--msp-color-text-muted);
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  flex: none;
}

.msp-viewer-volume-transfer-colorRow[aria-expanded="true"] .msp-viewer-volume-transfer-chevron {
  transform: rotate(180deg);
}

.msp-viewer-volume-transfer-palettes {
  min-width: 200px;
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-raised);
  box-shadow: inset 0 0 0 1px var(--msp-glass-edge),
    var(--msp-shadow-overlay);
  display: grid;
}

.msp-viewer-volume-transfer-palette {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 28px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  text-align: left;
  display: flex;

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    outline: none;
  }
}

.msp-viewer-volume-transfer-paletteSwatches {
  border-radius: var(--msp-radius-sm);
  width: 72px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
  flex: none;
  display: flex;
  overflow: hidden;
}

.msp-viewer-volume-transfer-paletteSwatch {
  flex: 1;
  height: 14px;
}

.msp-viewer-volume-transfer-paletteLabel {
  font-size: var(--msp-font-size-sm);
}
}
`);
