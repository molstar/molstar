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
  'actionFooter': 'msp-home-home-tab-actionFooter',
  'actionHeader': 'msp-home-home-tab-actionHeader',
  'actionTitle': 'msp-home-home-tab-actionTitle',
  'caption': 'msp-home-home-tab-caption',
  'captionRow': 'msp-home-home-tab-captionRow',
  'chevron': 'msp-home-home-tab-chevron',
  'drop': 'msp-home-home-tab-drop',
  'dropHint': 'msp-home-home-tab-dropHint',
  'dropText': 'msp-home-home-tab-dropText',
  'dropTitle': 'msp-home-home-tab-dropTitle',
  'group': 'msp-home-home-tab-group',
  'hint': 'msp-home-home-tab-hint',
  'home': 'msp-home-home-tab-home',
  'none': 'msp-home-home-tab-none',
  'recent': 'msp-home-home-tab-recent',
  'recentDescription': 'msp-home-home-tab-recentDescription',
  'recentLabel': 'msp-home-home-tab-recentLabel',
  'server': 'msp-home-home-tab-server',
  'spin': 'msp-home-home-tab-spin',
} as const;

defineStyles('views/home/home-tab.module.css', `@layer msp.components {
.msp-home-home-tab-home {
  align-content: start;
  gap: var(--msp-space-4);
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  display: grid;
}

.msp-home-home-tab-drop {
  align-items: center;
  gap: var(--msp-space-4);
  padding: var(--msp-space-4);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-color-control);
  color: var(--msp-color-text-muted);
  text-align: left;
  transition: background-color var(--msp-motion-fast) var(--msp-ease),
    color var(--msp-motion-fast) var(--msp-ease);
  display: flex;

  &:hover, &[data-dragging] {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[data-dragging] {
    outline: 2px dashed var(--msp-color-accent);
    outline-offset: -2px;
  }
}

.msp-home-home-tab-dropText {
  gap: 2px;
  min-width: 0;
  display: grid;
}

.msp-home-home-tab-dropTitle {
  color: var(--msp-color-text);
  font-weight: 600;
}

.msp-home-home-tab-dropHint {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
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

  &[data-open] {
    background: color-mix(in oklab, var(--msp-color-control), transparent 50%);
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

.msp-home-home-tab-actionFooter {
  justify-content: flex-end;
  display: flex;
}

.msp-home-home-tab-recent {
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

.msp-home-home-tab-recentLabel {
  color: var(--msp-color-text);
  font-weight: 500;
}

.msp-home-home-tab-recentDescription {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
}

.msp-home-home-tab-hint {
  align-items: center;
  gap: var(--msp-space-3);
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
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

.msp-home-home-tab-server {
  align-items: center;
  gap: var(--msp-space-3);
  padding: var(--msp-space-1) var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  grid-template-columns: auto minmax(0, 1fr);
  display: grid;
}
}
`);
