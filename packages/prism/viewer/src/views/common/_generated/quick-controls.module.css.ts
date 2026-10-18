/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from quick-controls.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'chevron': 'msp-viewer-quick-controls-chevron',
  'colors': 'msp-viewer-quick-controls-colors',
  'more': 'msp-viewer-quick-controls-more',
  'quick': 'msp-viewer-quick-controls-quick',
  'range': 'msp-viewer-quick-controls-range',
  'rangeLabel': 'msp-viewer-quick-controls-rangeLabel',
  'rangeValue': 'msp-viewer-quick-controls-rangeValue',
  'slider': 'msp-viewer-quick-controls-slider',
  'swatch': 'msp-viewer-quick-controls-swatch',
} as const;

defineStyles('viewer/views/common/quick-controls.module.css', `@layer msp.components {
.msp-viewer-quick-controls-colors {
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-quick-controls-quick {
  align-items: center;
  gap: var(--msp-space-1);
  display: flex;
}

.msp-viewer-quick-controls-swatch {
  border-radius: var(--msp-radius-sm);
  min-width: 0;
  max-width: 22px;
  height: 18px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  flex: 1 1 0;

  &:hover {
    transform: scale(1.08);
  }

  &[aria-selected="true"] {
    box-shadow: 0 0 0 1.5px var(--msp-color-bg),
      0 0 0 3px var(--msp-color-text);
  }
}

.msp-viewer-quick-controls-more {
  border-radius: var(--msp-radius-sm);
  width: 22px;
  height: 22px;
  color: var(--msp-color-text-muted);
  flex: none;
  place-items: center;
  margin-left: auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-viewer-quick-controls-chevron {
  transition: transform var(--msp-motion-fast) var(--msp-ease);
}

.msp-viewer-quick-controls-more[aria-expanded="true"] .msp-viewer-quick-controls-chevron {
  transform: rotate(180deg);
}

.msp-viewer-quick-controls-range {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 24px;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-viewer-quick-controls-rangeLabel {
  flex: none;
  min-width: 7ch;
}

.msp-viewer-quick-controls-slider {
  --msp-range-fill: color-mix(in oklab, var(--msp-color-text-muted), transparent 45%);
  -webkit-appearance: none;
  appearance: none;
  cursor: pointer;
  background: none;
  flex: 1;
  min-width: 0;
  height: 16px;

  &::-webkit-slider-runnable-track {
    background: linear-gradient(to right,
      var(--msp-range-fill) var(--msp-progress, 0%),
      var(--msp-color-pressed) var(--msp-progress, 0%));
    border-radius: 1.5px;
    height: 3px;
  }

  &::-webkit-slider-thumb {
    background: var(--msp-color-text-muted);
    -webkit-appearance: none;
    appearance: none;
    border-radius: 50%;
    width: 10px;
    height: 10px;
    margin-top: -3.5px;
  }

  &::-moz-range-track {
    background: var(--msp-color-pressed);
    border-radius: 1.5px;
    height: 3px;
  }

  &::-moz-range-progress {
    background: var(--msp-range-fill);
    border-radius: 1.5px;
    height: 3px;
  }

  &::-moz-range-thumb {
    background: var(--msp-color-text-muted);
    border: 0;
    border-radius: 50%;
    width: 10px;
    height: 10px;
  }
}

.msp-viewer-quick-controls-rangeValue {
  min-width: 4ch;
  color: var(--msp-color-text);
  font-variant-numeric: tabular-nums;
  text-align: right;
  flex: none;
}
}
`);
