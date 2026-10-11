/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from tooltip.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'muted': 'msp-tooltip-muted',
  'tooltip': 'msp-tooltip-tooltip',
  'tooltip-in': 'msp-tooltip-tooltip-in',
} as const;

defineStyles('components/tooltip/tooltip.module.css', `@layer msp.components {
.msp-tooltip-tooltip {
  z-index: var(--msp-z-tooltip);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-text);
  max-width: 240px;
  color: var(--msp-color-bg);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
  pointer-events: none;
  animation: msp-tooltip-tooltip-in var(--msp-motion-fast) var(--msp-ease);
  padding: 5px 9px;
  line-height: 1.3;
  position: fixed;
  top: 0;
  left: 0;

  &[hidden] {
    display: none;
  }
}

.msp-tooltip-tooltip[data-kind="row"] {
  border-radius: var(--msp-radius-lg);
  background: var(--msp-color-raised);
  max-width: 320px;
  color: var(--msp-color-text);
  box-shadow: 0 0 0 1px var(--msp-glass-edge),
    var(--msp-shadow-float);
  font-size: var(--msp-font-size-sm);
  white-space: normal;
  padding: 6px 10px;
  line-height: 1.4;
}

.msp-tooltip-muted {
  color: var(--msp-color-text-muted);
}

@keyframes msp-tooltip-tooltip-in {
  from {
    opacity: 0;
  }
}
}
`);
