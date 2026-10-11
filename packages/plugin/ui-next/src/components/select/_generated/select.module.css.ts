/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from select.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'chevron': 'msp-select-chevron',
  'native': 'msp-select-native',
  'select': 'msp-select-select',
} as const;

defineStyles('components/select/select.module.css', `@layer msp.components {
.msp-select-select {
  min-width: 0;
  height: var(--msp-control-height);
  display: inline-flex;
  position: relative;

  &[data-size="sm"] {
    height: 24px;
  }
}

.msp-select-native {
  min-width: 0;
  height: 100%;
  padding: 0 calc(var(--msp-space-3) + 18px) 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  -webkit-appearance: none;
  appearance: none;
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font: inherit;
  font-size: var(--msp-font-size-sm);
  text-overflow: ellipsis;
  cursor: pointer;
  transition: background-color var(--msp-motion-fast) var(--msp-ease);
  border: 0;
  flex: 1;
  margin: 0;

  &:hover {
    background: var(--msp-color-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 1px;
  }

  & option {
    background: var(--msp-color-panel);
    color: var(--msp-color-text);
  }
}

.msp-select-chevron {
  top: 50%;
  right: var(--msp-space-3);
  color: var(--msp-color-text-muted);
  pointer-events: none;
  position: absolute;
  transform: translateY(-50%);
}
}
`);
