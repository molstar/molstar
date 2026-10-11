/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from selection-bar.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'action': 'msp-canvas-selection-bar-action',
  'bar': 'msp-canvas-selection-bar-bar',
  'bar-in': 'msp-canvas-selection-bar-bar-in',
  'summary': 'msp-canvas-selection-bar-summary',
} as const;

defineStyles('views/canvas/selection-bar.module.css', `@layer msp.components {
.msp-canvas-selection-bar-bar {
  height: 34px;
  padding: 0 0 0 var(--msp-space-4);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  animation: msp-canvas-selection-bar-bar-in var(--msp-motion-base) var(--msp-ease);
  align-items: center;
  display: flex;
  overflow: hidden;
}

@keyframes msp-canvas-selection-bar-bar-in {
  from {
    opacity: 0;
    transform: translateY(4px);
  }
}

.msp-canvas-selection-bar-summary {
  align-items: center;
  gap: var(--msp-space-2);
  margin-right: var(--msp-space-3);
  font-size: var(--msp-font-size-sm);
  white-space: nowrap;
  display: flex;
}

.msp-canvas-selection-bar-action {
  width: 34px;
  height: 100%;
  color: var(--msp-color-text-muted);
  flex: none;
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

:where([data-docked]) .msp-canvas-selection-bar-bar {
  width: 100%;
  box-shadow: none;
  -webkit-backdrop-filter: none;
  border-top: 1px solid var(--msp-glass-edge);
  background: none;
  border-radius: 0;
  animation: none;
}

:where([data-docked]) .msp-canvas-selection-bar-summary {
  flex: 1;
}
}
`);
