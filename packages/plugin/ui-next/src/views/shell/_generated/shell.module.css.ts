/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from shell.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'canvas': 'msp-shell-canvas',
  'layer': 'msp-shell-layer',
  'shell': 'msp-shell-shell',
  'viewport': 'msp-shell-viewport',
} as const;

defineStyles('views/shell/shell.module.css', `@layer msp.components {
.msp-shell-shell {
  background: var(--msp-color-bg);
  grid-template: "top top top"
                 "left viewport right" minmax(0, 1fr)
                 "bottom bottom bottom"
                 / auto minmax(0, 1fr) auto;
  width: 100%;
  height: 100%;
  display: grid;
  position: relative;
  overflow: hidden;
  container-type: size;

  &[data-mode="compact"] {
    grid-template: "top"
                   "viewport" minmax(0, 1fr)
                   "bottom"
                   "right"
                   / minmax(0, 1fr);
  }
}

.msp-shell-viewport {
  grid-area: viewport;
  min-width: 0;
  min-height: 0;
  position: relative;
  overflow: hidden;
}

.msp-shell-canvas {
  position: absolute;
  inset: 0;
}

.msp-shell-layer {
  pointer-events: none;
  position: absolute;
  inset: 0;

  & > * {
    pointer-events: auto;
  }
}
}
`);
