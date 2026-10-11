/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from section-list.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'body': 'msp-shell-section-list-body',
  'chevron': 'msp-shell-section-list-chevron',
  'header': 'msp-shell-section-list-header',
  'list': 'msp-shell-section-list-list',
  'section': 'msp-shell-section-list-section',
  'title': 'msp-shell-section-list-title',
} as const;

defineStyles('views/shell/section-list.module.css', `@layer msp.components {
.msp-shell-section-list-list {
  gap: var(--msp-space-1);
  padding: 0 var(--msp-space-3) var(--msp-space-3);
  display: grid;
}

.msp-shell-section-list-section {
  border-radius: var(--msp-radius-md);

  &:has([aria-expanded="true"]) {
    background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  }
}

.msp-shell-section-list-header {
  align-items: center;
  gap: var(--msp-space-3);
  width: 100%;
  height: var(--msp-row-height);
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  text-align: left;
  font-weight: 600;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-shell-section-list-title {
  flex: 1;
}

.msp-shell-section-list-chevron {
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="false"] > & {
    transform: rotate(-90deg);
  }
}

.msp-shell-section-list-body {
  padding: 0 var(--msp-space-3) var(--msp-space-3);
}

.msp-shell-section-list-list[data-flat] {
  padding: 0;

  & .msp-shell-section-list-section {
    background: none;
  }

  & .msp-shell-section-list-header {
    height: 26px;
    font-size: var(--msp-font-size-sm);
    padding: 0;

    &:hover {
      color: var(--msp-color-accent);
      background: none;
    }
  }

  & .msp-shell-section-list-body {
    padding: 0 0 var(--msp-space-2);
  }
}
}
`);
