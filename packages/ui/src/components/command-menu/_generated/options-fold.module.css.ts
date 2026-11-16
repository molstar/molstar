/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from options-fold.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'body': 'msui-command-menu-options-fold-body',
  'chevron': 'msui-command-menu-options-fold-chevron',
  'fold': 'msui-command-menu-options-fold-fold',
  'head': 'msui-command-menu-options-fold-head',
  'summary': 'msui-command-menu-options-fold-summary',
  'title': 'msui-command-menu-options-fold-title',
} as const;

defineStyles('components/command-menu/options-fold.module.css', `@layer msui.components {
.msui-command-menu-options-fold-fold {
  border-top: 1px solid var(--msui-line);
  padding: var(--msui-space-1) var(--msui-surface-pad);
}

.msui-command-menu-options-fold-head {
  grid-template-columns: var(--msui-icon-slot) auto minmax(0, 1fr);
  column-gap: var(--msui-icon-gap);
  width: 100%;
  height: var(--msui-row);
  padding-inline: var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text);
  text-align: left;
  align-items: center;
  display: grid;

  &:hover {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-command-menu-options-fold-chevron {
  width: var(--msui-icon-slot);
  color: var(--msui-text-muted);
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  justify-content: center;
  align-items: center;
  display: inline-flex;

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}

.msui-command-menu-options-fold-title {
  white-space: nowrap;
}

.msui-command-menu-options-fold-summary {
  color: var(--msui-text-muted);
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  padding-inline-start: var(--msui-space-2);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-command-menu-options-fold-body {
  padding-block: var(--msui-space-1);
}
}
`);
