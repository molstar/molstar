/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from structure-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'actions': 'msp-structure-structure-tab-actions',
  'fact': 'msp-structure-structure-tab-fact',
  'facts': 'msp-structure-structure-tab-facts',
  'inspector': 'msp-structure-structure-tab-inspector',
  'inspectorHeader': 'msp-structure-structure-tab-inspectorHeader',
  'inspectorKind': 'msp-structure-structure-tab-inspectorKind',
  'inspectorTitle': 'msp-structure-structure-tab-inspectorTitle',
  'tab': 'msp-structure-structure-tab-tab',
} as const;

defineStyles('views/structure/structure-tab.module.css', `@layer msp.components {
.msp-structure-structure-tab-tab {
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-4);
  font-size: var(--msp-font-size-sm);
  display: grid;
}

.msp-structure-structure-tab-actions {
  gap: var(--msp-space-1);
  display: flex;
}

.msp-structure-structure-tab-inspector {
  gap: var(--msp-space-2);
  padding: var(--msp-space-2) var(--msp-space-3) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 50%);
  display: grid;
}

.msp-structure-structure-tab-inspectorHeader {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 26px;
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

.msp-structure-structure-tab-facts {
  gap: var(--msp-space-1) var(--msp-space-4);
  font-size: var(--msp-font-size-xs);
  flex-wrap: wrap;
  margin: 0;
  display: flex;
}

.msp-structure-structure-tab-fact {
  gap: var(--msp-space-2);
  display: flex;

  & dt {
    color: var(--msp-color-text-subtle);
  }

  & dd {
    font-variant-numeric: tabular-nums;
    margin: 0;
  }
}
}
`);
