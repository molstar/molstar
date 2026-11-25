/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from button.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'block': 'msui-button-block',
  'button': 'msui-button-button',
  'danger': 'msui-button-danger',
  'ghost': 'msui-button-ghost',
  'icon': 'msui-button-icon',
  'iconOnly': 'msui-button-iconOnly',
  'primary': 'msui-button-primary',
  'small': 'msui-button-small',
  'tiny': 'msui-button-tiny',
} as const;

defineStyles('components/button/button.module.css', `@layer msui.components {
.msui-button-button {
  --msui-button-height: var(--msui-row);
  --msui-button-icon: var(--msui-icon);
  justify-content: center;
  align-items: center;
  gap: var(--msui-space-1);
  height: max(var(--msui-button-height), var(--msui-touch-target, 0px));
  padding: 0 var(--msui-space-2);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  white-space: nowrap;
  transition: background var(--msui-motion-fast) var(--msui-ease),
    color var(--msui-motion-fast) var(--msui-ease);
  display: inline-flex;

  &:hover:not(:disabled) {
    background: var(--msui-control-hover);
  }

  &:active:not(:disabled), &[aria-pressed="true"] {
    background: var(--msui-control-active);
    box-shadow: inset 0 0 0 1px var(--msui-selected-edge);
  }

  &:disabled {
    color: var(--msui-text-subtle);
  }
}

.msui-button-primary {
  background: color-mix(in oklab, var(--msui-accent), transparent 84%);
  color: var(--msui-accent-strong);
  font-weight: var(--msui-weight-strong);

  &:hover:not(:disabled) {
    background: color-mix(in oklab, var(--msui-accent), transparent 76%);
  }

  &:active:not(:disabled), &[aria-pressed="true"] {
    background: color-mix(in oklab, var(--msui-accent), transparent 70%);
    box-shadow: inset 0 0 0 1px var(--msui-selected-edge);
  }

  &:disabled {
    background: color-mix(in oklab, var(--msui-accent), transparent 92%);
    color: var(--msui-text-subtle);
  }
}

.msui-button-block {
  width: 100%;
}

.msui-button-ghost {
  background: none;
}

.msui-button-danger:hover:not(:disabled) {
  background: var(--msui-danger);
  color: var(--msui-on-accent);
}

.msui-button-small {
  --msui-button-height: var(--msui-control-height);
  --msui-button-icon: var(--msui-icon-sm);
  padding: 0 var(--msui-control-inset);
}

.msui-button-tiny {
  --msui-button-height: var(--msui-control-height);
  --msui-button-icon: var(--msui-icon-sm);
  padding: 0 var(--msui-space-1);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
}

.msui-button-iconOnly {
  width: max(var(--msui-button-height), var(--msui-touch-target, 0px));
  padding: 0;
}

.msui-button-icon {
  --msui-icon-size: var(--msui-button-icon);
  display: inline-flex;
}
}
`);
