/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from hover-label.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'label': 'msp-canvas-hover-label-label',
  'primary': 'msp-canvas-hover-label-primary',
  'secondary': 'msp-canvas-hover-label-secondary',
} as const;

defineStyles('views/canvas/hover-label.module.css', `@layer msp.components {
.msp-canvas-hover-label-label {
  z-index: var(--msp-z-tooltip);
  max-width: 280px;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  font-size: var(--msp-font-size-sm);
  pointer-events: none;
  will-change: transform;
  gap: 1px;
  display: grid;
  position: absolute;
  top: 0;
  left: 0;

  &[hidden] {
    display: none;
  }

  background: var(--msp-glass-bg);
  border: var(--msp-glass-border);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
}

.msp-canvas-hover-label-primary {
  font-weight: 600;
}

.msp-canvas-hover-label-secondary {
  color: var(--msp-color-text-muted);
}
}
`);
