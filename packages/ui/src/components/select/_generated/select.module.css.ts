/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from select.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'button': 'msui-select-button',
  'check': 'msui-select-check',
  'chevron': 'msui-select-chevron',
  'floating': 'msui-select-floating',
  'group': 'msui-select-group',
  'groupBody': 'msui-select-groupBody',
  'groupChevron': 'msui-select-groupChevron',
  'groupHead': 'msui-select-groupHead',
  'label': 'msui-select-label',
  'list': 'msui-select-list',
  'listContent': 'msui-select-listContent',
  'menuPopover': 'msui-select-menuPopover',
  'option': 'msui-select-option',
  'optionLabel': 'msui-select-optionLabel',
  'select': 'msui-select-select',
} as const;

defineStyles('components/select/select.module.css', `@layer msui.components {
.msui-select-select {
  gap: var(--msui-space-1);
  min-width: 0;
  display: grid;
}

.msui-select-button {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: center;
  gap: var(--msui-icon-gap);
  width: 100%;
  min-width: 0;
  height: var(--msui-row);
  padding: 0 var(--msui-space-1) 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  text-align: left;
  display: flex;

  &[data-size="sm"] {
    height: var(--msui-control-height);
  }

  &:hover, &[aria-expanded="true"] {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-select-label {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-select-chevron {
  color: var(--msui-text-muted);
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  flex: none;

  [aria-expanded="true"] > & {
    transform: rotate(180deg);
  }
}

.msui-select-list {
  --msui-icon-size: var(--msui-icon-sm);
  border-radius: var(--msui-radius-lg);
  background: color-mix(in oklab, var(--msui-control), transparent 40%);
  max-height: 240px;
  box-shadow: inset 0 0 0 1px var(--msui-line);
}

.msui-select-listContent {
  padding: var(--msui-space-1) var(--msui-surface-pad);
  display: grid;
}

.msui-select-option {
  justify-content: space-between;
  align-items: center;
  gap: var(--msui-space-2);
  min-height: var(--msui-row);
  padding-inline: calc(var(--msui-control-inset) - var(--msui-surface-pad))
    calc(var(--msui-space-1) - var(--msui-surface-pad));
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  text-align: left;
  display: flex;

  &:hover, &:focus-visible {
    background: var(--msui-control-hover);
    color: var(--msui-text);
    outline: none;
  }

  &[aria-selected="true"] {
    color: var(--msui-accent-strong);
  }
}

.msui-select-optionLabel {
  min-width: 0;
}

.msui-select-group {
  display: grid;
}

.msui-select-groupHead {
  justify-content: space-between;
  align-items: center;
  gap: var(--msui-space-2);
  min-height: var(--msui-row);
  padding-inline: calc(var(--msui-control-inset) - var(--msui-surface-pad))
    calc(var(--msui-space-1) - var(--msui-surface-pad));
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  font-weight: var(--msui-weight-strong);
  text-align: left;
  display: flex;

  &:hover, &:focus-visible {
    background: var(--msui-control-hover);
    color: var(--msui-text);
    outline: none;
  }
}

.msui-select-groupChevron {
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  flex: none;

  [aria-expanded="true"] > & {
    transform: rotate(180deg);
  }
}

.msui-select-groupBody {
  display: grid;
}

.msui-select-check {
  visibility: hidden;
  flex: none;

  [aria-selected="true"] > & {
    visibility: visible;
  }
}

.msui-select-floating {
  color: inherit;
  width: max-content;
  max-width: calc(100vw - 16px);
  animation: msui-float-in var(--msui-enter);
  background: none;
  border: 0;
  margin: 0;
  padding: 0;
  position: fixed;
  inset: auto;
  overflow: visible;

  &[data-side="below"] {
    --msui-rise-y: calc(-1 * var(--msui-rise));
  }

  &[data-side="above"] {
    --msui-rise-y: var(--msui-rise);
  }

  &[data-side="right"] {
    --msui-rise-x: calc(-1 * var(--msui-rise));
  }

  &[data-side="left"] {
    --msui-rise-x: var(--msui-rise);
  }

  & .msui-select-option {
    white-space: nowrap;
  }

  & > .msui-select-list {
    background: var(--msui-surface-raised);
    box-shadow: inset 0 0 0 1px var(--msui-line),
      var(--msui-shadow-overlay);
  }
}

.msui-select-menuPopover {
  --msui-icon-size: var(--msui-icon-sm);
  min-width: 180px;
  max-height: calc(100vh - 16px);
  padding: var(--msui-space-1) var(--msui-surface-pad);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-overlay);
  font-size: var(--msui-text-sm);
  display: grid;
  overflow-y: auto;
}
}
`);
