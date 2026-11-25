/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from settings-view.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'search': 'msui-forms-settings-view-search',
  'settings': 'msui-forms-settings-view-settings',
  'title': 'msui-forms-settings-view-title',
} as const;

defineStyles('forms/settings-view.module.css', `@layer msui.components {
.msui-forms-settings-view-settings {
  min-width: 0;
}

.msui-forms-settings-view-search {
  flex: auto;
  min-width: 0;
}

.msui-forms-settings-view-title {
  --msui-icon-size: var(--msui-icon);
  align-items: center;
  gap: var(--msui-space-2);
  color: var(--msui-text);
  display: inline-flex;

  & svg {
    color: var(--msui-text-muted);
    flex: none;
  }
}
}
`);
