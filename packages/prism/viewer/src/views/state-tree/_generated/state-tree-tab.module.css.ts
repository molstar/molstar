/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from state-tree-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'blocks': 'msp-viewer-state-tree-tab-blocks',
  'caption': 'msp-viewer-state-tree-tab-caption',
  'error': 'msp-viewer-state-tree-tab-error',
  'transform': 'msp-viewer-state-tree-tab-transform',
  'tree': 'msp-viewer-state-tree-tab-tree',
} as const;

defineStyles('viewer/views/state-tree/state-tree-tab.module.css', `@layer msp.components {
.msp-viewer-state-tree-tab-error {
  color: var(--msp-color-danger);
}

.msp-viewer-state-tree-tab-transform {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
}

.msp-viewer-state-tree-tab-blocks {
  gap: var(--msp-space-1);
  margin: 0 calc(-1 * var(--msp-space-2));
  display: grid;
}

.msp-viewer-state-tree-tab-caption {
  margin: var(--msp-space-2) 0 0;
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-weight: 500;
}

.msp-viewer-state-tree-tab-tree {
  --msp-tree-indent: 8px;
}
}
`);
