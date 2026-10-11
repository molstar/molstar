/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from palette.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'before': 'msp-palette-before',
  'caption': 'msp-palette-caption',
  'command': 'msp-palette-command',
  'empty': 'msp-palette-empty',
  'error': 'msp-palette-error',
  'footer': 'msp-palette-footer',
  'head': 'msp-palette-head',
  'headIcon': 'msp-palette-headIcon',
  'hint': 'msp-palette-hint',
  'input': 'msp-palette-input',
  'keys': 'msp-palette-keys',
  'list': 'msp-palette-list',
  'message': 'msp-palette-message',
  'more': 'msp-palette-more',
  'palette': 'msp-palette-palette',
  'palette-fade': 'msp-palette-palette-fade',
  'palette-in': 'msp-palette-palette-in',
  'row': 'msp-palette-row',
  'rowIcon': 'msp-palette-rowIcon',
  'rowMeta': 'msp-palette-rowMeta',
  'rowTitle': 'msp-palette-rowTitle',
  'run': 'msp-palette-run',
  'scrim': 'msp-palette-scrim',
  'sentence': 'msp-palette-sentence',
  'spacer': 'msp-palette-spacer',
  'tail': 'msp-palette-tail',
  'token': 'msp-palette-token',
  'value': 'msp-palette-value',
} as const;

defineStyles('views/palette/palette.module.css', `@layer msp.components {
.msp-palette-scrim {
  z-index: var(--msp-z-dialog);
  padding: 12vh var(--msp-space-4) var(--msp-space-4);
  background: color-mix(in oklab, var(--msp-color-bg), transparent 75%);
  animation: msp-palette-palette-fade var(--msp-motion-fast) var(--msp-ease);
  justify-content: center;
  align-items: flex-start;
  display: flex;
  position: absolute;
  inset: 0;
}

@keyframes msp-palette-palette-fade {
  from {
    opacity: 0;
  }
}

.msp-palette-palette {
  border-radius: var(--msp-radius-lg);
  background: var(--msp-color-panel);
  width: min(600px, 100%);
  max-height: min(520px, 76vh);
  box-shadow: var(--msp-glass-shadow),
    0 12px 32px #0000001f;
  font-size: var(--msp-font-size-md);
  animation: msp-palette-palette-in var(--msp-motion-base) var(--msp-ease);
  flex-direction: column;
  display: flex;
  overflow: hidden;
}

@keyframes msp-palette-palette-in {
  from {
    opacity: 0;
    transform: translateY(-6px);
  }
}

.msp-palette-head {
  align-items: center;
  gap: var(--msp-space-3);
  height: 48px;
  padding: 0 var(--msp-space-4);
  border-bottom: 1px solid var(--msp-glass-edge);
  flex: none;
  display: flex;
}

.msp-palette-headIcon {
  color: var(--msp-color-text-subtle);
  flex: none;
}

.msp-palette-input {
  min-width: 0;
  height: 100%;
  color: var(--msp-color-text);
  font: inherit;
  background: none;
  border: 0;
  outline: none;
  flex: 1;

  &::placeholder {
    color: var(--msp-color-text-subtle);
  }
}

.msp-palette-command {
  align-items: center;
  gap: var(--msp-space-2);
  height: 26px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-accent), transparent 85%);
  color: var(--msp-color-accent);
  font-size: var(--msp-font-size-sm);
  flex: none;
  font-weight: 600;
  display: flex;
}

.msp-palette-message, .msp-palette-error {
  padding: var(--msp-space-2) var(--msp-space-4);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  margin: 0;
}

.msp-palette-error {
  color: var(--msp-color-danger);
}

.msp-palette-hint {
  color: var(--msp-color-text-muted);
}

.msp-palette-list {
  min-height: 0;
  padding: var(--msp-space-2);
  flex: 1;
  overflow: auto;

  &:empty {
    display: none;
  }
}

.msp-palette-empty {
  padding: var(--msp-space-5) var(--msp-space-4);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-sm);
  text-align: center;
  margin: 0;
}

.msp-palette-caption {
  padding: var(--msp-space-3) var(--msp-space-3) var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-palette-row {
  align-items: center;
  gap: var(--msp-space-3);
  height: 34px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text);
  cursor: default;
  grid-template-columns: 20px minmax(0, 1fr) auto auto;
  display: grid;

  &[aria-selected="true"] {
    background: var(--msp-color-hover);
  }

  &[aria-disabled="true"] {
    color: var(--msp-color-text-subtle);
  }
}

.msp-palette-rowIcon {
  color: var(--msp-color-text-muted);
  place-items: center;
  display: grid;
}

.msp-palette-rowTitle {
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;

  & mark {
    color: var(--msp-color-accent);
    background: none;
    font-weight: 600;
  }
}

.msp-palette-tail {
  margin-left: var(--msp-space-2);
  color: var(--msp-color-text-muted);
}

.msp-palette-rowMeta {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
}

.msp-palette-keys {
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  color: var(--msp-color-text-muted);
  font: inherit;
  font-size: var(--msp-font-size-xs);
  padding: 1px 6px;
}

.msp-palette-sentence {
  align-items: center;
  gap: var(--msp-space-2);
  padding: var(--msp-space-3) var(--msp-space-4);
  border-bottom: 1px solid var(--msp-glass-edge);
  flex-wrap: wrap;
  flex: none;
  display: flex;
}

.msp-palette-token {
  align-items: center;
  gap: var(--msp-space-2);
  display: inline-flex;
}

.msp-palette-before {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-sm);
}

.msp-palette-value {
  height: 26px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);

  &:hover {
    background: var(--msp-color-hover);
  }

  &[data-empty] {
    color: var(--msp-color-text-subtle);
  }

  &[data-missing] {
    box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-danger), transparent 50%);
  }

  &[data-active] {
    box-shadow: inset 0 0 0 1.5px var(--msp-color-accent);
  }
}

.msp-palette-spacer {
  flex: 1;
}

.msp-palette-run {
  align-items: center;
  gap: var(--msp-space-2);
  height: 28px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-accent);
  color: #fff;
  font-size: var(--msp-font-size-sm);
  font-weight: 600;
  display: flex;

  & kbd {
    font: inherit;
    opacity: .8;
  }
}

.msp-palette-more {
  padding: var(--msp-space-2) var(--msp-space-4);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-palette-footer {
  gap: var(--msp-space-4);
  padding: var(--msp-space-2) var(--msp-space-4);
  border-top: 1px solid var(--msp-glass-edge);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
  display: flex;

  & kbd {
    border-radius: var(--msp-radius-sm);
    background: var(--msp-color-control);
    font: inherit;
    padding: 0 4px;
  }
}
}
`);
