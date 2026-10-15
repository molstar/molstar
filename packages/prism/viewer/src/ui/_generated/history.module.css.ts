/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from history.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'bar': 'msp-viewer-history-bar',
  'count': 'msp-viewer-history-count',
  'detail': 'msp-viewer-history-detail',
  'empty': 'msp-viewer-history-empty',
  'history': 'msp-viewer-history-history',
  'label': 'msp-viewer-history-label',
  'origin': 'msp-viewer-history-origin',
  'step': 'msp-viewer-history-step',
  'time': 'msp-viewer-history-time',
  'unsafe': 'msp-viewer-history-unsafe',
} as const;

defineStyles('viewer/ui/history.module.css', `@layer msp.components {
.msp-viewer-history-history {
  gap: var(--msp-space-2);
  padding: 0 var(--msp-panel-inset) var(--msp-panel-inset);
  font-size: var(--msp-font-size-sm);
  display: grid;
}

.msp-viewer-history-bar {
  align-items: center;
  gap: var(--msp-space-3);
  display: flex;
}

.msp-viewer-history-count {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin-left: auto;
}

.msp-viewer-history-empty {
  padding: var(--msp-space-4);
  color: var(--msp-color-text-muted);
}

.msp-viewer-history-step {
  align-items: baseline;
  gap: var(--msp-space-1) var(--msp-space-3);
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  grid-template-columns: auto minmax(0, 1fr) auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[data-unsafe] .msp-viewer-history-label {
    color: var(--msp-color-text-muted);
  }
}

.msp-viewer-history-origin {
  color: var(--msp-color-text-subtle);
  align-self: center;
}

.msp-viewer-history-label {
  overflow-wrap: anywhere;
}

.msp-viewer-history-time {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
}

.msp-viewer-history-detail {
  align-items: center;
  gap: var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  grid-column: 2 / -1;
  display: flex;
}

.msp-viewer-history-unsafe {
  overflow-wrap: anywhere;
  opacity: .8;
  min-width: 0;
}
}
`);
