/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from popover.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'popover': 'msui-popover-popover',
  'title': 'msui-popover-title',
} as const;

defineStyles('components/popover/popover.module.css', `@layer msui.components {
.msui-popover-popover {
  width: max-content;
  max-width: min(340px, 100vw - 16px);
  max-height: calc(100vh - 16px);
  padding: var(--msui-space-1) var(--msui-surface-pad) var(--msui-space-2);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-overlay);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  align-content: start;
  display: grid;
  overflow-y: auto;

  [data-column] > & {
    width: auto;
    max-width: none;
  }
}

.msui-popover-title {
  height: var(--msui-header);
  padding-inline: var(--msui-inner);
  font-weight: var(--msui-weight-strong);
  color: var(--msui-text);
  align-items: center;
  display: flex;
}
}
`);
