/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from superpose.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'entries': 'msp-structure-superpose-entries',
  'entry': 'msp-structure-superpose-entry',
  'hint': 'msp-structure-superpose-hint',
  'index': 'msp-structure-superpose-index',
  'label': 'msp-structure-superpose-label',
  'role': 'msp-structure-superpose-role',
  'superpose': 'msp-structure-superpose-superpose',
} as const;

defineStyles('views/structure/superpose.module.css', `@layer msp.components {
.msp-structure-superpose-superpose {
  font-size: var(--msp-font-size-sm);
  gap: var(--msp-space-2);
  display: grid;
}

.msp-structure-superpose-entries {
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  display: grid;
}

.msp-structure-superpose-hint {
  padding: var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  line-height: 1.45;
}

.msp-structure-superpose-entry {
  align-items: center;
  gap: var(--msp-space-2);
  height: 26px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  text-align: left;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-structure-superpose-index {
  min-width: 1.5ch;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  flex: none;
}

.msp-structure-superpose-label {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  font-weight: 500;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-structure-superpose-role {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
}
}
`);
