/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from toolbar.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'button': 'msui-toolbar-button',
  'caption': 'msui-toolbar-caption',
  'check': 'msui-toolbar-check',
  'chevron': 'msui-toolbar-chevron',
  'divider': 'msui-toolbar-divider',
  'end': 'msui-toolbar-end',
  'face': 'msui-toolbar-face',
  'group': 'msui-toolbar-group',
  'heading': 'msui-toolbar-heading',
  'hint': 'msui-toolbar-hint',
  'label': 'msui-toolbar-label',
  'level': 'msui-toolbar-level',
  'mark': 'msui-toolbar-mark',
  'menu': 'msui-toolbar-menu',
  'menuRoot': 'msui-toolbar-menuRoot',
  'row': 'msui-toolbar-row',
  'slot': 'msui-toolbar-slot',
  'toolbar': 'msui-toolbar-toolbar',
} as const;

defineStyles('components/toolbar/toolbar.module.css', `@layer msui.components {
.msui-toolbar-toolbar {
  align-items: flex-start;
  gap: var(--msui-space-2);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  flex-direction: row;
  display: flex;

  &[data-orientation="vertical"] {
    flex-direction: column;
    align-items: stretch;
    width: max-content;
  }
}

.msui-toolbar-group {
  --msui-icon-size: var(--msui-icon);
  padding: var(--msui-surface-pad);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface);
  box-shadow: var(--msui-surface-edge);
  flex-direction: row;
  align-items: center;
  display: flex;

  [data-orientation="vertical"] > & {
    flex-direction: column;
    align-items: stretch;
  }

  [data-elevated] > & {
    background: var(--msui-surface-raised);
    box-shadow: inset 0 0 0 1px var(--msui-line),
      var(--msui-shadow-float);
  }

  &:not(:has( > .msui-toolbar-button:not([hidden]))) {
    display: none;
  }
}

.msui-toolbar-toolbar[data-grouping="separators"] {
  padding: var(--msui-surface-pad);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface);
  box-shadow: var(--msui-surface-edge);
  gap: 0;

  &[data-elevated] {
    background: var(--msui-surface-raised);
    box-shadow: inset 0 0 0 1px var(--msui-line),
      var(--msui-shadow-float);
  }

  & > .msui-toolbar-group {
    box-shadow: none;
    background: none;
    border-radius: 0;
    padding: 0;
  }

  & > .msui-toolbar-group:has( > .msui-toolbar-button:not([hidden])) ~ .msui-toolbar-group:before {
    content: "";
    width: 1px;
    height: var(--msui-icon);
    margin-inline: var(--msui-space-1);
    background: var(--msui-line-strong);
    flex: none;
    align-self: center;
  }

  &[data-orientation="vertical"] > .msui-toolbar-group:has( > .msui-toolbar-button:not([hidden])) ~ .msui-toolbar-group:before {
    width: var(--msui-icon);
    height: 1px;
    margin: var(--msui-space-1) 0;
  }
}

.msui-toolbar-button {
  justify-content: center;
  align-items: center;
  gap: var(--msui-space-1);
  min-width: var(--msui-header);
  height: var(--msui-header);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  white-space: nowrap;
  transition: background var(--msui-motion-fast) var(--msui-ease),
    color var(--msui-motion-fast) var(--msui-ease);
  flex: none;
  padding: 0;
  display: inline-flex;
  position: relative;

  &:hover:not([aria-disabled="true"]), &[aria-expanded="true"] {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &:active:not([aria-disabled="true"]) {
    background: var(--msui-control-active);
  }

  &[aria-pressed="true"] {
    background: var(--msui-control-active);
    box-shadow: inset 0 0 0 1px var(--msui-selected-edge);
    color: var(--msui-accent-strong);
  }

  &[aria-disabled="true"] {
    color: var(--msui-text-subtle);
    cursor: default;
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &[hidden] {
    display: none;
  }
}

[data-captions="beside"] .msui-toolbar-button[data-captioned] {
  padding-inline: var(--msui-space-2);
}

[data-captions="below"] .msui-toolbar-button[data-captioned] {
  height: calc(var(--msui-header) + var(--msui-space-2));
  padding: var(--msui-space-1);
  flex-direction: column;
  gap: 0;
}

.msui-toolbar-face {
  display: inline-flex;
}

.msui-toolbar-caption {
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  font-variant-numeric: tabular-nums;
}

.msui-toolbar-mark {
  right: var(--msui-space-1);
  bottom: var(--msui-space-1);
  border-top: var(--msui-space-1) solid transparent;
  border-right: var(--msui-space-1) solid currentColor;
  opacity: .5;
  pointer-events: none;
  position: absolute;
}

.msui-toolbar-menuRoot, .msui-toolbar-level {
  display: grid;
}

.msui-toolbar-menu {
  --msui-icon-size: var(--msui-icon-sm);
  min-width: calc(40 * var(--msui-unit));
  max-height: calc(100vh - var(--msui-space-4));
  padding: var(--msui-space-1) var(--msui-surface-pad);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-overlay);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  display: grid;
  overflow-y: auto;
}

.msui-toolbar-row {
  grid-template-columns: var(--msui-icon-slot) 1fr auto;
  align-items: center;
  column-gap: var(--msui-icon-gap);
  height: var(--msui-row);
  padding-inline: var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  text-align: left;
  white-space: nowrap;
  display: grid;

  [data-ends] > & {
    grid-template-columns: var(--msui-icon-slot) 1fr auto var(--msui-icon-sm);
  }

  &:hover:not([aria-disabled="true"]), &:focus-visible, &[aria-expanded="true"] {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &[role="menuitemradio"][aria-checked="true"] {
    color: var(--msui-accent-strong);
  }

  &[aria-disabled="true"] {
    color: var(--msui-text-subtle);
    cursor: default;
  }
}

.msui-toolbar-slot {
  width: var(--msui-icon-slot);
  justify-content: center;
  align-items: center;
  display: inline-flex;
}

.msui-toolbar-label {
  white-space: nowrap;
}

.msui-toolbar-hint {
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  padding-inline-start: var(--msui-space-2);

  &:empty {
    padding: 0;
  }
}

.msui-toolbar-end {
  justify-content: center;
  align-items: center;
  display: inline-flex;
}

.msui-toolbar-check {
  visibility: hidden;

  [aria-checked="true"] > .msui-toolbar-end > & {
    visibility: visible;
  }
}

.msui-toolbar-chevron {
  color: var(--msui-text-subtle);
}

.msui-toolbar-divider {
  height: 1px;
  margin: var(--msui-space-1) var(--msui-inner);
  background: var(--msui-line);
}

.msui-toolbar-heading {
  padding: var(--msui-space-2) var(--msui-inner) var(--msui-space-1) calc(var(--msui-inner) + var(--msui-indent));
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);

  &[data-first] {
    padding-top: var(--msui-space-1);
  }
}
}
`);
