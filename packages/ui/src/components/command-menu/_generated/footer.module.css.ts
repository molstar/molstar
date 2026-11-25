/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from footer.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'cap': 'msui-command-menu-footer-cap',
  'footer': 'msui-command-menu-footer-footer',
  'hint': 'msui-command-menu-footer-hint',
  'keys': 'msui-command-menu-footer-keys',
} as const;

defineStyles('components/command-menu/footer.module.css', `@layer msui.components {
.msui-command-menu-footer-footer {
  align-items: center;
  column-gap: var(--msui-space-4);
  row-gap: var(--msui-space-1);
  min-height: var(--msui-header);
  padding: var(--msui-space-1) var(--msui-inner);
  box-shadow: 0 -1px 0 var(--msui-line);
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  flex-wrap: wrap;
  flex: none;
  display: flex;
  position: relative;
}

.msui-command-menu-footer-hint {
  align-items: center;
  gap: var(--msui-space-1);
  white-space: nowrap;
  display: inline-flex;
}

.msui-command-menu-footer-keys {
  font: inherit;
  gap: 1px;
  display: inline-flex;
}

.msui-command-menu-footer-cap {
  min-width: var(--msui-icon);
  height: var(--msui-icon);
  padding-inline: var(--msui-space-1);
  border-radius: var(--msui-radius-sm);
  background: var(--msui-surface-sunken);
  box-shadow: inset 0 0 0 1px var(--msui-line);
  color: var(--msui-text-muted);
  font-family: var(--msui-font-sans);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  white-space: nowrap;
  place-items: center;
  display: inline-grid;
}
}
`);
