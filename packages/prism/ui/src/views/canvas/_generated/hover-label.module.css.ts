/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from hover-label.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'code': 'msp-canvas-hover-label-code',
  'corner': 'msp-canvas-hover-label-corner',
  'item': 'msp-canvas-hover-label-item',
  'label': 'msp-canvas-hover-label-label',
  'link': 'msp-canvas-hover-label-link',
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

  background: var(--msp-glass-panel-bg);
  box-shadow: var(--msp-glass-shadow);
}

.msp-canvas-hover-label-corner {
  top: auto;
  left: auto;
  right: var(--msp-canvas-right-inset, var(--msp-canvas-edge, var(--msp-space-3)));
  bottom: calc(var(--msp-canvas-edge, var(--msp-space-3)) + var(--msp-corner-room, 0px));
  will-change: auto;
}

.msp-canvas-hover-label-primary {
  font-weight: 600;
}

.msp-canvas-hover-label-secondary {
  color: var(--msp-color-text-muted);
}

.msp-canvas-hover-label-item {
  text-indent: -1em;
  padding-inline-start: 1em;
}

.msp-canvas-hover-label-link {
  color: var(--msp-color-accent);
  text-decoration: underline;
}

.msp-canvas-hover-label-code {
  font-family: var(--msp-font-mono, ui-monospace, monospace);
  font-size: .95em;
}
}
`);
