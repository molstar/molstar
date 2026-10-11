/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from home-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'action': 'msp-home-home-tab-action',
  'actionBody': 'msp-home-home-tab-actionBody',
  'actionDescription': 'msp-home-home-tab-actionDescription',
  'actionHeader': 'msp-home-home-tab-actionHeader',
  'actionTitle': 'msp-home-home-tab-actionTitle',
  'caption': 'msp-home-home-tab-caption',
  'captionActions': 'msp-home-home-tab-captionActions',
  'captionRow': 'msp-home-home-tab-captionRow',
  'chevron': 'msp-home-home-tab-chevron',
  'content': 'msp-home-home-tab-content',
  'group': 'msp-home-home-tab-group',
  'hint': 'msp-home-home-tab-hint',
  'home': 'msp-home-home-tab-home',
  'none': 'msp-home-home-tab-none',
  'server': 'msp-home-home-tab-server',
  'spin': 'msp-home-home-tab-spin',
  'state': 'msp-home-home-tab-state',
  'stateDate': 'msp-home-home-tab-stateDate',
  'stateLabel': 'msp-home-home-tab-stateLabel',
  'states': 'msp-home-home-tab-states',
} as const;

defineStyles('views/home/home-tab.module.css', `@layer msp.components {
.msp-home-home-tab-home {
  height: 100%;
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  display: flex;
}

.msp-home-home-tab-content {
  align-content: start;
  gap: var(--msp-space-4);
  min-height: 0;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  flex: 1 1 0;
  display: grid;
  overflow: auto;
}

.msp-home-home-tab-states {
  grid-template-columns: minmax(0, 1fr);
  max-height: 264px;
  display: grid;
  overflow: auto;
}

.msp-home-home-tab-group {
  gap: var(--msp-space-1);
  display: grid;
}

.msp-home-home-tab-caption {
  padding: 0 var(--msp-space-2) var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  font-weight: 500;
}

.msp-home-home-tab-action {
  border-radius: var(--msp-radius-md);
  overflow: hidden;

  &[data-open] {
    background: var(--msp-color-well);
  }
}

.msp-home-home-tab-actionHeader {
  align-items: center;
  gap: var(--msp-space-3);
  width: 100%;
  height: 32px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  text-align: left;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-home-home-tab-actionTitle {
  color: var(--msp-color-text);
  flex: none;
  font-weight: 600;
}

.msp-home-home-tab-actionDescription {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  overflow: hidden;
}

.msp-home-home-tab-chevron {
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="false"] > & {
    transform: rotate(-90deg);
  }
}

.msp-home-home-tab-actionBody {
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) var(--msp-space-3) var(--msp-space-3);
  display: grid;
}

.msp-home-home-tab-state {
  align-items: center;
  gap: var(--msp-space-3);
  height: 30px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  text-align: left;
  grid-template-columns: auto auto minmax(0, 1fr);
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-home-home-tab-stateLabel {
  color: var(--msp-color-text);
  font-weight: 500;
}

.msp-home-home-tab-stateDate {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
}

.msp-home-home-tab-hint {
  margin: 0 var(--msp-space-2) var(--msp-space-2);
  align-items: center;
  gap: var(--msp-space-3);
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  & kbd {
    border-radius: var(--msp-radius-sm);
    background: var(--msp-color-control);
    color: var(--msp-color-text-muted);
    font: inherit;
    padding: 1px 6px;
  }
}

.msp-home-home-tab-spin {
  animation: .9s linear infinite msp-home-home-tab-spin;
}

@keyframes msp-home-home-tab-spin {
  to {
    transform: rotate(360deg);
  }
}

.msp-home-home-tab-captionRow {
  padding-right: var(--msp-space-1);
  justify-content: space-between;
  align-items: center;
  display: flex;
}

.msp-home-home-tab-none {
  padding: var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-home-home-tab-captionActions {
  gap: 2px;
  display: flex;
}

.msp-home-home-tab-server {
  align-items: center;
  gap: var(--msp-space-3);
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  grid-template-columns: auto minmax(0, 1fr);
  display: grid;
}
}
`);
