/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from log-view.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'entries': 'msui-log-log-view-entries',
  'filters': 'msui-log-log-view-filters',
  'kind': 'msui-log-log-view-kind',
  'kinds': 'msui-log-log-view-kinds',
  'log': 'msui-log-log-view-log',
  'row': 'msui-log-log-view-row',
  'search': 'msui-log-log-view-search',
} as const;

defineStyles('log/log-view.module.css', `@layer msui.components {
.msui-log-log-view-log {
  gap: var(--msui-space-1);
  flex-direction: column;
  flex: 1 1 0;
  min-width: 0;
  min-height: 0;
  display: flex;
}

.msui-log-log-view-filters {
  flex: none;
}

.msui-log-log-view-kinds {
  flex: none;
  display: flex;
}

.msui-log-log-view-search {
  flex: auto;
  min-width: 0;
}

.msui-log-log-view-entries {
  flex: 1 1 0;
  min-height: 0;
}

.msui-log-log-view-row {
  &[data-kind="error"] .msui-log-log-view-kind {
    color: var(--msui-danger);
  }

  &[data-kind="warning"] .msui-log-log-view-kind {
    color: var(--msui-warning);
  }

  &[data-kind="success"] .msui-log-log-view-kind {
    color: var(--msui-success);
  }

  &[data-kind="debug"] {
    color: var(--msui-text-muted);
  }
}

.msui-log-log-view-kind {
  color: var(--msui-text-subtle);
}
}
`);
