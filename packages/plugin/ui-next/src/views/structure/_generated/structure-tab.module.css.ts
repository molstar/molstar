/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from structure-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'actions': 'msp-structure-structure-tab-actions',
  'inspector': 'msp-structure-structure-tab-inspector',
  'inspectorHeader': 'msp-structure-structure-tab-inspectorHeader',
  'inspectorKind': 'msp-structure-structure-tab-inspectorKind',
  'inspectorTitle': 'msp-structure-structure-tab-inspectorTitle',
  'tab': 'msp-structure-structure-tab-tab',
} as const;

defineStyles('views/structure/structure-tab.module.css', `@layer msp.components {
.msp-structure-structure-tab-tab {
  gap: var(--msp-space-3);
  padding: var(--msp-space-2) var(--msp-space-2) var(--msp-space-4);
  display: grid;
}

.msp-structure-structure-tab-actions {
  gap: var(--msp-space-1);
  display: flex;
}

.msp-structure-structure-tab-inspector {
  gap: var(--msp-space-2);
  padding: var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 50%);
  display: grid;
}

.msp-structure-structure-tab-inspectorHeader {
  align-items: center;
  gap: var(--msp-space-2);
  padding-bottom: var(--msp-space-1);
  display: flex;
}

.msp-structure-structure-tab-inspectorTitle {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  min-width: 0;
  font-weight: 600;
  overflow: hidden;
}

.msp-structure-structure-tab-inspectorKind {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-transform: capitalize;
}
}
`);
