/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from text-input.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'input': 'msui-text-input-input',
} as const;

defineStyles('components/text-input/text-input.module.css', `@layer msui.components {
.msui-text-input-input {
  min-width: 0;
  height: var(--msui-row);
  padding: 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  color: var(--msui-text);
  font: inherit;
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  border: 0;
  margin: 0;

  &[data-size="sm"] {
    height: var(--msui-control-height);
  }

  &:hover:not(:disabled) {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &::placeholder, &:disabled {
    color: var(--msui-text-subtle);
  }
}
}
`);
