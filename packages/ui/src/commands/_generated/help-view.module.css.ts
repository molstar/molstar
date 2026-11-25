/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from help-view.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'action': 'msui-commands-help-view-action',
  'conflict': 'msui-commands-help-view-conflict',
  'conflicts': 'msui-commands-help-view-conflicts',
  'edit': 'msui-commands-help-view-edit',
  'editBar': 'msui-commands-help-view-editBar',
  'group': 'msui-commands-help-view-group',
  'headingRow': 'msui-commands-help-view-headingRow',
  'help': 'msui-commands-help-view-help',
  'input': 'msui-commands-help-view-input',
  'key': 'msui-commands-help-view-key',
  'keys': 'msui-commands-help-view-keys',
  'list': 'msui-commands-help-view-list',
  'or': 'msui-commands-help-view-or',
  'recorder': 'msui-commands-help-view-recorder',
  'row': 'msui-commands-help-view-row',
  'srOnly': 'msui-commands-help-view-srOnly',
} as const;

defineStyles('commands/help-view.module.css', `@layer msui.components {
.msui-commands-help-view-help {
  min-width: 0;
}

.msui-commands-help-view-group {
  gap: var(--msui-space-1);
  min-width: 0;
  display: grid;
}

.msui-commands-help-view-list {
  margin: 0;
  display: grid;
}

.msui-commands-help-view-row {
  height: auto;
  min-height: var(--msui-row);
  column-gap: var(--msui-space-3);

  &:hover {
    background: var(--msui-control-hover);
  }
}

.msui-commands-help-view-action {
  min-width: 0;
}

.msui-commands-help-view-input {
  justify-content: flex-end;
  align-items: center;
  gap: var(--msui-space-1) var(--msui-space-2);
  flex-wrap: wrap;
  margin: 0;
  display: flex;
}

.msui-commands-help-view-keys {
  gap: var(--msui-space-1);
  font: inherit;
  display: inline-flex;
}

.msui-commands-help-view-key {
  min-width: var(--msui-icon);
  padding: 0 var(--msui-space-1);
  border-radius: var(--msui-radius-sm);
  background: var(--msui-control);
  color: var(--msui-text);
  font: inherit;
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  text-align: center;
  white-space: nowrap;
}

.msui-commands-help-view-or {
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
}

.msui-commands-help-view-conflicts {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: flex-start;
  gap: var(--msui-space-2);
  padding: var(--msui-space-2) var(--msui-inner);
  border-radius: var(--msui-radius-md);
  background: color-mix(in oklab, var(--msui-warning), transparent 85%);
  color: var(--msui-text);
  display: flex;

  & > svg {
    margin-top: var(--msui-space-1);
    color: var(--msui-warning);
    flex: none;
  }
}

.msui-commands-help-view-conflict {
  margin: 0;
}

.msui-commands-help-view-headingRow {
  justify-content: space-between;
  align-items: flex-end;
  gap: var(--msui-space-2);
  padding-inline-end: var(--msui-inner);
  display: flex;
}

.msui-commands-help-view-edit {
  gap: var(--msui-space-1);
  display: inline-flex;
}

.msui-commands-help-view-recorder {
  padding: 0 var(--msui-space-2);
  border-radius: var(--msui-radius-sm);
  outline: 1px solid var(--msui-accent);
  color: var(--msui-accent);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  white-space: nowrap;
}

.msui-commands-help-view-editBar {
  justify-content: space-between;
  align-items: center;
  gap: var(--msui-space-2);
  padding: var(--msui-space-1) var(--msui-inner) 0;
  display: flex;
}

.msui-commands-help-view-srOnly {
  clip-path: inset(50%);
  white-space: nowrap;
  width: 1px;
  height: 1px;
  position: absolute;
  overflow: hidden;
}
}
`);
