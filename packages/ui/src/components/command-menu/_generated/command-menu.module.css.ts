/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from command-menu.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'ask': 'msui-command-menu-ask',
  'askKeys': 'msui-command-menu-askKeys',
  'askLabel': 'msui-command-menu-askLabel',
  'askSlot': 'msui-command-menu-askSlot',
  'askText': 'msui-command-menu-askText',
  'caption': 'msui-command-menu-caption',
  'crumb': 'msui-command-menu-crumb',
  'empty': 'msui-command-menu-empty',
  'head': 'msui-command-menu-head',
  'headLine': 'msui-command-menu-headLine',
  'headSlot': 'msui-command-menu-headSlot',
  'input': 'msui-command-menu-input',
  'list': 'msui-command-menu-list',
  'listContent': 'msui-command-menu-listContent',
  'menu': 'msui-command-menu-menu',
  'message': 'msui-command-menu-message',
  'rest': 'msui-command-menu-rest',
  'row': 'msui-command-menu-row',
  'rowEnd': 'msui-command-menu-rowEnd',
  'rowMeta': 'msui-command-menu-rowMeta',
  'rowSlot': 'msui-command-menu-rowSlot',
  'rowTitle': 'msui-command-menu-rowTitle',
  'scrim': 'msui-command-menu-scrim',
  'shared': 'msui-command-menu-shared',
  'sharedOption': 'msui-command-menu-sharedOption',
  'sharedTitle': 'msui-command-menu-sharedTitle',
  'srOnly': 'msui-command-menu-srOnly',
  'tail': 'msui-command-menu-tail',
} as const;

defineStyles('components/command-menu/command-menu.module.css', `@layer msui.components {
.msui-command-menu-scrim {
  z-index: var(--msui-z-dialog);
  padding: 12vh var(--msui-space-4) var(--msui-space-4);
  background: color-mix(in oklab, var(--msui-bg), transparent 70%);
  justify-content: center;
  align-items: flex-start;
  display: flex;
  position: absolute;
  inset: 0;
}

.msui-command-menu-menu {
  --msui-icon-size: var(--msui-icon-sm);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface);
  width: min(560px, 100%);
  max-height: min(480px, 76vh);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-overlay);
  color: var(--msui-text);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
  flex-direction: column;
  display: flex;
  overflow: hidden;
}

.msui-command-menu-head {
  grid-template-columns: var(--msui-icon-slot) minmax(0, 1fr);
  column-gap: var(--msui-icon-gap);
  min-height: calc(10 * var(--msui-unit));
  padding: var(--msui-space-2) var(--msui-inset);
  border-bottom: 1px solid var(--msui-line);
  flex: none;
  align-items: start;
  display: grid;
}

.msui-command-menu-headSlot {
  --msui-icon-size: var(--msui-icon);
  width: var(--msui-icon-slot);
  height: var(--msui-row);
  color: var(--msui-text-muted);
  justify-content: center;
  align-items: center;
  display: inline-flex;
}

.msui-command-menu-headLine {
  align-items: center;
  gap: var(--msui-space-1);
  min-width: 0;
  height: var(--msui-row);
  display: flex;
}

.msui-command-menu-input {
  min-width: 0;
  height: var(--msui-row);
  color: var(--msui-text);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
  background: none;
  border: 0;
  outline: none;
  flex: 1;
  padding: 0;

  &::placeholder {
    color: var(--msui-text-subtle);
  }
}

.msui-command-menu-crumb {
  align-items: center;
  gap: var(--msui-space-1);
  height: var(--msui-control-height);
  padding-inline: var(--msui-space-1);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  white-space: nowrap;
  flex: none;
  display: inline-flex;

  &:hover {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }
}

.msui-command-menu-message {
  padding: var(--msui-space-1) var(--msui-inset);
  color: var(--msui-text-muted);
  margin: 0;
}

.msui-command-menu-list {
  flex: auto;
  min-height: 0;

  &:has(.msui-command-menu-listContent:empty) {
    display: none;
  }
}

.msui-command-menu-listContent {
  padding: var(--msui-space-1) var(--msui-surface-pad);
  align-content: start;
  display: grid;
}

.msui-command-menu-empty {
  padding: var(--msui-space-5) var(--msui-inset);
  color: var(--msui-text-subtle);
  text-align: center;
  margin: 0;
}

.msui-command-menu-rest {
  pointer-events: none;
}

.msui-command-menu-caption {
  height: var(--msui-row);
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  align-items: flex-end;
  padding-block-end: var(--msui-space-1);
  padding-inline-start: calc(var(--msui-inner) + var(--msui-icon-slot) + var(--msui-icon-gap));
  display: flex;
}

.msui-command-menu-row {
  grid-template-columns: var(--msui-icon-slot) minmax(max-content, 1fr) minmax(0, auto) auto;
  column-gap: var(--msui-icon-gap);
  height: var(--msui-row);
  padding-inline: var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text);
  cursor: default;
  align-items: center;
  display: grid;

  &[aria-selected="true"] {
    background: var(--msui-control-hover);
  }

  &[aria-selected="true"] .msui-command-menu-rowMeta {
    color: var(--msui-text-muted);
  }

  &[aria-disabled="true"] {
    color: var(--msui-text-subtle);
  }
}

.msui-command-menu-rowSlot {
  width: var(--msui-icon-slot);
  height: var(--msui-icon);
  color: var(--msui-text-muted);
  justify-content: center;
  align-items: center;
  display: inline-flex;
}

.msui-command-menu-rowTitle {
  white-space: nowrap;

  & mark {
    color: var(--msui-accent);
    font-weight: var(--msui-weight-strong);
    background: none;
  }

  &[data-current] {
    color: var(--msui-accent);
  }
}

.msui-command-menu-tail {
  color: var(--msui-text-muted);
  margin-inline-start: var(--msui-space-1);
}

.msui-command-menu-rowMeta {
  text-align: end;
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  justify-self: end;
  padding-inline-start: var(--msui-space-2);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-command-menu-rowEnd {
  color: var(--msui-text-subtle);
  align-items: center;
  margin-inline-start: var(--msui-space-2);
  display: inline-flex;
}

.msui-command-menu-shared {
  align-items: center;
  gap: var(--msui-space-1);
  min-height: var(--msui-header);
  padding-inline: calc(var(--msui-inset) + var(--msui-icon-slot) + var(--msui-icon-gap)) var(--msui-inset);
  border-bottom: 1px solid var(--msui-line);
  flex: none;
  display: flex;
}

.msui-command-menu-sharedTitle {
  color: var(--msui-text-subtle);
  margin-inline-end: var(--msui-space-2);
}

.msui-command-menu-sharedOption {
  height: var(--msui-control-height);
  padding-inline: var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  white-space: nowrap;

  &:hover {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &[aria-pressed="true"] {
    background: var(--msui-control);
    color: var(--msui-accent);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-command-menu-ask {
  grid-template-columns: var(--msui-icon-slot) auto minmax(0, 1fr) auto;
  column-gap: var(--msui-icon-gap);
  height: var(--msui-row);
  margin: 0 var(--msui-surface-pad) var(--msui-space-1);
  padding-inline: var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text);
  text-align: left;
  flex: none;
  align-items: center;
  display: grid;

  &:hover, &[data-default] {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &[data-setup] {
    grid-template-columns: var(--msui-icon-slot) minmax(0, 1fr);
    color: var(--msui-text-muted);
  }
}

.msui-command-menu-askSlot {
  width: var(--msui-icon-slot);
  color: var(--msui-accent);
  justify-content: center;
  align-items: center;
  display: inline-flex;

  [data-setup] > & {
    color: var(--msui-text-subtle);
  }
}

.msui-command-menu-askLabel {
  white-space: nowrap;
}

.msui-command-menu-askText {
  color: var(--msui-text-muted);
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  padding-inline-start: var(--msui-space-2);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-command-menu-askKeys {
  color: var(--msui-text-subtle);
  font: inherit;
  font-size: var(--msui-text-xs);
}

.msui-command-menu-srOnly {
  clip-path: inset(50%);
  white-space: nowrap;
  width: 1px;
  height: 1px;
  position: absolute;
  overflow: hidden;
}
}
`);
