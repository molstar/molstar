/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from context-menu.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'check': 'msui-context-menu-check',
  'keys': 'msui-context-menu-keys',
  'label': 'msui-context-menu-label',
  'menu': 'msui-context-menu-menu',
  'row': 'msui-context-menu-row',
  'slot': 'msui-context-menu-slot',
} as const;

defineStyles('components/context-menu/context-menu.module.css', `@layer msui.components {
.msui-context-menu-menu {
  --msui-icon-size: var(--msui-icon-sm);
  z-index: var(--msui-z-overlay);
  min-width: 180px;
  padding: var(--msui-space-1) var(--msui-surface-pad);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-overlay);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  --msui-rise-x: 0px;
  --msui-rise-y: calc(-1 * var(--msui-rise));
  animation: msui-float-in var(--msui-motion-fast) var(--msui-ease);
  display: grid;
  position: fixed;
  top: 0;
  left: 0;
}

.msui-context-menu-row {
  align-items: center;
  column-gap: var(--msui-icon-gap);
  height: var(--msui-row);
  padding-inline: var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  text-align: left;
  grid-template-columns: 1fr auto;
  display: grid;

  [data-icons] > & {
    grid-template-columns: var(--msui-icon-slot) 1fr auto;
  }

  &:hover:not([aria-disabled="true"]), &[aria-expanded="true"], &:focus-visible {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &[role="menuitemradio"][aria-checked="true"] {
    color: var(--msui-accent-strong);
  }

  &[aria-disabled="true"] {
    color: var(--msui-text-subtle);
    cursor: default;
  }
}

.msui-context-menu-row[data-separated] {
  margin-top: calc(var(--msui-space-1) * 2 + 1px);
  position: relative;

  &:before {
    content: "";
    inset-inline: calc(-1 * var(--msui-surface-pad));
    top: calc(-1 * var(--msui-space-1) - 1px);
    background: var(--msui-line);
    height: 1px;
    position: absolute;
  }
}

.msui-context-menu-slot {
  width: var(--msui-icon-slot);
  height: var(--msui-icon);
  justify-content: center;
  align-items: center;
  display: inline-flex;
}

.msui-context-menu-label {
  white-space: nowrap;
}

.msui-context-menu-keys {
  color: var(--msui-text-subtle);
  font-family: inherit;
  font-size: var(--msui-text-xs);
  padding-inline-start: var(--msui-space-2);
}

.msui-context-menu-check {
  --msui-icon-size: var(--msui-icon-sm);
  width: var(--msui-icon-slot);
  justify-content: center;
  align-items: center;
  padding-inline-start: var(--msui-space-2);
  display: inline-flex;
}
}
`);
