/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from history-view.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'action': 'msui-history-history-view-action',
  'blocked': 'msui-history-history-view-blocked',
  'count': 'msui-history-history-view-count',
  'history': 'msui-history-history-view-history',
  'row': 'msui-history-history-view-row',
  'steps': 'msui-history-history-view-steps',
  'unsafe': 'msui-history-history-view-unsafe',
} as const;

defineStyles('history/history-view.module.css', `@layer msui.components {
.msui-history-history-view-history {
  gap: var(--msui-space-1);
  flex-direction: column;
  flex: 1 1 0;
  min-width: 0;
  min-height: 0;
  display: flex;
}

.msui-history-history-view-count {
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  font-variant-numeric: tabular-nums;
  text-align: end;
  flex: 1;
}

.msui-history-history-view-blocked {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: center;
  gap: var(--msui-space-1);
  min-height: var(--msui-row);
  padding-inline: var(--msui-inner);
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  margin: 0;
  display: flex;
}

.msui-history-history-view-steps {
  flex: 1 1 0;
  min-height: 0;
}

.msui-history-history-view-row {
  &[data-state="current"] {
    background: var(--msui-accent-muted);

    & .msui-slot {
      color: var(--msui-accent);
    }
  }

  &[data-state="undone"] {
    color: var(--msui-text-subtle);
  }

  &[data-state="run"] {
    color: var(--msui-text-muted);
  }

  &:hover {
    background: var(--msui-control-hover);
  }

  &:not(:hover, :focus-within) .msui-history-history-view-action {
    opacity: 0;
  }
}

.msui-history-history-view-action {
  flex: none;
}

.msui-history-history-view-unsafe {
  --msui-icon-size: var(--msui-icon-sm);
  color: var(--msui-warning);
  display: inline-flex;
}
}
`);
