/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from choices.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'choice': 'msui-choices-choice',
  'choices': 'msui-choices-choices',
  'content': 'msui-choices-content',
  'label': 'msui-choices-label',
  'row': 'msui-choices-row',
} as const;

defineStyles('components/choices/choices.module.css', `@layer msui.components {
.msui-choices-choices {
  grid-template-columns: minmax(0, 1fr);
  min-width: 0;
  display: grid;
}

.msui-choices-choices[data-labelled] {
  grid-template-columns: var(--msui-label) minmax(0, 1fr);
  column-gap: var(--msui-space-3);
}

.msui-form > .msui-choices-choices {
  grid-column: 1 / -1;
  grid-template-columns: subgrid;

  &:not([data-labelled]) > .msui-choices-row {
    grid-column: 2;
  }
}

.msui-choices-label {
  min-width: 0;
  height: var(--msui-row);
  color: var(--msui-text-muted);
  white-space: nowrap;
  align-self: start;
  align-items: center;
  display: flex;
  overflow: hidden;
}

.msui-choices-row {
  --msui-choices-min: 64px;
  grid-template-columns: repeat(auto-fit, minmax(min(var(--msui-choices-min), 100%), 1fr));
  grid-auto-rows: var(--msui-row);
  align-items: center;
  column-gap: var(--msui-space-1);
  min-width: 0;
  display: grid;

  &[data-icons] {
    --msui-choices-min: calc(var(--msui-control-height) + 2 * var(--msui-space-1));
  }

  &[data-columns] {
    grid-template-columns: repeat(var(--msui-choices-columns), minmax(0, 1fr));
  }
}

.msui-choices-choice {
  --msui-icon-size: var(--msui-icon-sm);
  min-width: 0;
  height: var(--msui-control-height);
  padding: 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  white-space: nowrap;
  justify-content: center;
  align-items: center;
  display: flex;

  &:hover:not([aria-disabled="true"]) {
    background: var(--msui-control-hover);
  }

  &:active:not([aria-disabled="true"]) {
    background: var(--msui-control-active);
  }

  &[aria-checked="true"] {
    background: color-mix(in oklab, var(--msui-accent), transparent 84%);
    color: var(--msui-accent-strong);
    box-shadow: inset 0 0 0 1px var(--msui-selected-edge);

    &:hover {
      background: color-mix(in oklab, var(--msui-accent), transparent 76%);
    }
  }

  &[aria-disabled="true"] {
    color: var(--msui-text-subtle);
    cursor: default;
  }

  [data-icons] > & {
    padding: 0;
  }
}

.msui-choices-content {
  align-items: center;
  gap: var(--msui-icon-gap);
  flex: none;
  display: inline-flex;
}
}
`);
