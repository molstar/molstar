/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from popover.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'popover': 'msp-popover-popover',
  'title': 'msp-popover-title',
} as const;

defineStyles('components/popover/popover.module.css', `@layer msp.components {
.msp-popover-popover {
  gap: var(--msp-space-2);
  width: min(340px, 100vw - 16px);
  max-height: calc(100vh - 16px);
  padding: var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-raised);
  box-shadow: inset 0 0 0 1px var(--msp-glass-edge),
    var(--msp-shadow-overlay);
  font-size: var(--msp-font-size-sm);
  display: grid;
  overflow-y: auto;
}

.msp-popover-title {
  color: var(--msp-color-text);
  font-weight: 600;
}
}
`);
