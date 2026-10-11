/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from button.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'button': 'msp-button-button',
  'danger': 'msp-button-danger',
  'ghost': 'msp-button-ghost',
  'icon': 'msp-button-icon',
  'iconOnly': 'msp-button-iconOnly',
  'primary': 'msp-button-primary',
  'small': 'msp-button-small',
} as const;

defineStyles('components/button/button.module.css', `@layer msp.components {
.msp-button-button {
  justify-content: center;
  align-items: center;
  gap: var(--msp-space-2);
  height: var(--msp-control-height);
  padding: 0 var(--msp-space-4);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-md);
  white-space: nowrap;
  transition: background-color var(--msp-motion-fast) var(--msp-ease),
    color var(--msp-motion-fast) var(--msp-ease);
  display: inline-flex;

  &:hover:not(:disabled) {
    background: var(--msp-color-hover);
  }

  &:active:not(:disabled), &[aria-pressed="true"] {
    background: var(--msp-color-pressed);
  }

  &:disabled {
    color: var(--msp-color-text-subtle);
  }
}

.msp-button-primary {
  color: var(--msp-color-accent);
  font-weight: 600;

  &:hover:not(:disabled) {
    color: var(--msp-color-accent-strong);
  }
}

.msp-button-ghost {
  background: none;
}

.msp-button-danger:hover:not(:disabled) {
  background: var(--msp-color-danger);
  color: #fff;
}

.msp-button-small {
  height: 24px;
  padding: 0 var(--msp-space-3);
  font-size: var(--msp-font-size-sm);
}

.msp-button-iconOnly {
  width: var(--msp-control-height);
  padding: 0;

  &.msp-button-small {
    width: 24px;
  }
}

.msp-button-icon {
  width: 18px;
  height: 18px;
  display: inline-flex;

  & > svg {
    width: 100%;
    height: 100%;
  }
}
}
`);
