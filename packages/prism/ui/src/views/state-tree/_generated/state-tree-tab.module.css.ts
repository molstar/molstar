/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from state-tree-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'blocks': 'msp-state-tree-state-tree-tab-blocks',
  'caption': 'msp-state-tree-state-tree-tab-caption',
  'error': 'msp-state-tree-state-tree-tab-error',
  'transform': 'msp-state-tree-state-tree-tab-transform',
  'tree': 'msp-state-tree-state-tree-tab-tree',
} as const;

defineStyles('views/state-tree/state-tree-tab.module.css', `@layer msp.components {
.msp-state-tree-state-tree-tab-error {
  color: var(--msp-color-danger);
}

.msp-state-tree-state-tree-tab-transform {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
}

.msp-state-tree-state-tree-tab-blocks {
  gap: var(--msp-space-1);
  margin: 0 calc(-1 * var(--msp-space-2));
  display: grid;
}

.msp-state-tree-state-tree-tab-caption {
  margin: var(--msp-space-2) 0 0;
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-weight: 500;
}

.msp-state-tree-state-tree-tab-tree {
  --msp-tree-indent: 8px;
}
}
`);
