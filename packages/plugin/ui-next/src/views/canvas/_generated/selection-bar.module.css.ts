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
  'dot': 'msp-canvas-selection-bar-dot',
  'summary': 'msp-canvas-selection-bar-summary',
} as const;

defineStyles('views/canvas/selection-bar.module.css', `@layer msp.components {
.msp-canvas-selection-bar-bar {
  height: 34px;
  padding: 0 3px 0 var(--msp-space-4);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  animation: msp-canvas-selection-bar-bar-in var(--msp-motion-base) var(--msp-ease);
  align-items: center;
  gap: 2px;
  display: flex;
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

.msp-canvas-selection-bar-dot {
  background: var(--msp-color-accent);
  border-radius: 50%;
  width: 8px;
  height: 8px;
}

.msp-canvas-selection-bar-action {
  border-radius: var(--msp-radius-md);
  width: 28px;
  height: 28px;
  color: var(--msp-color-text-muted);
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}
}
`);
