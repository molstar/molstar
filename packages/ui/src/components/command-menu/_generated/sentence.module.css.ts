/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sentence.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'before': 'msui-command-menu-sentence-before',
  'below': 'msui-command-menu-sentence-below',
  'command': 'msui-command-menu-sentence-command',
  'edit': 'msui-command-menu-sentence-edit',
  'end': 'msui-command-menu-sentence-end',
  'field': 'msui-command-menu-sentence-field',
  'fieldView': 'msui-command-menu-sentence-fieldView',
  'form': 'msui-command-menu-sentence-form',
  'input': 'msui-command-menu-sentence-input',
  'line': 'msui-command-menu-sentence-line',
  'lineHint': 'msui-command-menu-sentence-lineHint',
  'proposed': 'msui-command-menu-sentence-proposed',
  'range': 'msui-command-menu-sentence-range',
  'rangeEnd': 'msui-command-menu-sentence-rangeEnd',
  'rangeRow': 'msui-command-menu-sentence-rangeRow',
  'rangeValue': 'msui-command-menu-sentence-rangeValue',
  'run': 'msui-command-menu-sentence-run',
  'runKeys': 'msui-command-menu-sentence-runKeys',
  'sentence': 'msui-command-menu-sentence-sentence',
  'slot': 'msui-command-menu-sentence-slot',
  'tail': 'msui-command-menu-sentence-tail',
  'token': 'msui-command-menu-sentence-token',
  'value': 'msui-command-menu-sentence-value',
} as const;

defineStyles('components/command-menu/sentence.module.css', `@layer msui.components {
.msui-command-menu-sentence-sentence {
  align-items: center;
  gap: var(--msui-space-1);
  min-width: 0;
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  flex-wrap: wrap;
  display: flex;
}

.msui-command-menu-sentence-command {
  height: var(--msui-row);
  font-weight: var(--msui-weight-strong);
  white-space: nowrap;
  align-items: center;
  margin-inline-end: var(--msui-space-1);
  display: inline-flex;
}

.msui-command-menu-sentence-field {
  align-items: center;
  gap: var(--msui-space-1);
  display: inline-flex;
}

.msui-command-menu-sentence-before {
  color: var(--msui-text-subtle);
  white-space: nowrap;
}

.msui-command-menu-sentence-token {
  height: var(--msui-row);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  align-items: center;
  display: inline-flex;

  &[data-state="missing"] {
    box-shadow: none;
    outline: 1px dashed var(--msui-line-strong);
    outline-offset: -1px;
    color: var(--msui-text-subtle);
    background: none;
    font-style: italic;
  }

  &[data-state="default"] {
    color: var(--msui-text-subtle);
  }

  &[data-proposed] {
    background: var(--msui-accent-muted);
  }

  &[data-active] {
    background: var(--msui-surface-sunken);
    box-shadow: inset 0 0 0 1px var(--msui-focus);
    outline: 0;
    font-style: normal;
  }

  &[data-error] {
    box-shadow: inset 0 0 0 1px var(--msui-danger);
  }
}

.msui-command-menu-sentence-value {
  height: 100%;
  padding-inline: var(--msui-control-inset);
  border-radius: inherit;
  color: inherit;
  font: inherit;
  white-space: nowrap;

  &:hover {
    background: var(--msui-control-hover);
  }
}

.msui-command-menu-sentence-slot, .msui-command-menu-sentence-tail {
  display: contents;
}

.msui-command-menu-sentence-edit {
  height: var(--msui-row);
  padding-inline: var(--msui-control-inset);
  align-items: center;
  display: inline-grid;

  &:after {
    content: attr(data-text) " ";
    visibility: hidden;
    white-space: pre;
    grid-area: 1 / 1;
  }

  .msui-command-menu-sentence-tail > & {
    min-width: calc(30 * var(--msui-unit));
    flex: 1;
    padding-inline: 0;
  }
}

.msui-command-menu-sentence-input {
  width: 100%;
  min-width: var(--msui-space-1);
  color: var(--msui-text);
  font: inherit;
  background: none;
  border: 0;
  outline: none;
  grid-area: 1 / 1;
  padding: 0;

  &::placeholder {
    color: var(--msui-text-subtle);
  }

  [data-state="set"] > &::placeholder {
    color: var(--msui-text);
  }

  [data-state="missing"] > &::placeholder {
    font-style: italic;
  }
}

.msui-command-menu-sentence-end {
  align-items: center;
  gap: var(--msui-space-1);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  margin-inline-start: auto;
  display: inline-flex;
}

.msui-command-menu-sentence-run {
  gap: var(--msui-space-2);
}

.msui-command-menu-sentence-runKeys > kbd {
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msui-on-accent), transparent 50%);
  color: inherit;
  background: none;
}

.msui-command-menu-sentence-line {
  padding: var(--msui-space-1) var(--msui-inset) 0
    calc(var(--msui-inset) + var(--msui-icon-slot) + var(--msui-icon-gap));
  color: var(--msui-text-muted);
  flex: none;
  margin: 0;

  &[data-tone="error"] {
    color: var(--msui-danger);
  }
}

.msui-command-menu-sentence-lineHint {
  color: var(--msui-text-muted);
}

.msui-command-menu-sentence-rangeRow {
  grid-template-columns: auto minmax(0, 1fr) auto calc(10 * var(--msui-unit));
  column-gap: var(--msui-space-2);
  height: var(--msui-row);
  padding-inline: calc(var(--msui-inset) + var(--msui-icon-slot) + var(--msui-icon-gap)) var(--msui-inset);
  flex: none;
  align-items: center;
  margin-block-start: var(--msui-space-1);
  display: grid;
}

.msui-command-menu-sentence-range {
  width: 100%;
  accent-color: var(--msui-accent);
  margin: 0;
}

.msui-command-menu-sentence-rangeEnd {
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
}

.msui-command-menu-sentence-rangeValue {
  color: var(--msui-text);
}

.msui-command-menu-sentence-fieldView {
  padding: var(--msui-space-1) var(--msui-inset);
  flex: none;
}

.msui-command-menu-sentence-below {
  min-height: 0;
  box-shadow: 0 -1px 0 var(--msui-line);
  flex: 0 auto;
  position: relative;
  overflow-y: auto;

  & > :first-child {
    border-top: 0;
  }

  &:empty {
    display: none;
  }
}

.msui-command-menu-sentence-form {
  border-top: 1px solid var(--msui-line);
}

.msui-command-menu-sentence-proposed {
  align-items: center;
  gap: var(--msui-icon-gap);
  padding: var(--msui-space-1) var(--msui-inset);
  background: var(--msui-accent-muted);
  color: var(--msui-text-muted);
  flex: none;
  margin: 0;
  display: flex;
}
}
`);
