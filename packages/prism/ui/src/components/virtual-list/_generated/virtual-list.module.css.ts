/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from virtual-list.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'content': 'msp-virtual-list-content',
  'list': 'msp-virtual-list-list',
  'window': 'msp-virtual-list-window',
} as const;

defineStyles('components/virtual-list/virtual-list.module.css', `@layer msp.components {
.msp-virtual-list-list {
  flex: 1 1 0;
}

.msp-virtual-list-content {
  contain: strict;
  width: 100%;
  position: relative;
}

.msp-virtual-list-window {
  will-change: transform;
  position: absolute;
  inset: 0 0 auto;
}
}
`);
