/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from console.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'action': 'msp-command-language-console-action',
  'bar': 'msp-command-language-console-bar',
  'card': 'msp-command-language-console-card',
  'chip': 'msp-command-language-console-chip',
  'code': 'msp-command-language-console-code',
  'command': 'msp-command-language-console-command',
  'completion': 'msp-command-language-console-completion',
  'completions': 'msp-command-language-console-completions',
  'composer': 'msp-command-language-console-composer',
  'console': 'msp-command-language-console-console',
  'detail': 'msp-command-language-console-detail',
  'does': 'msp-command-language-console-does',
  'drop': 'msp-command-language-console-drop',
  'empty': 'msp-command-language-console-empty',
  'emptyText': 'msp-command-language-console-emptyText',
  'emptyTitle': 'msp-command-language-console-emptyTitle',
  'entry': 'msp-command-language-console-entry',
  'example': 'msp-command-language-console-example',
  'examples': 'msp-command-language-console-examples',
  'field': 'msp-command-language-console-field',
  'hint': 'msp-command-language-console-hint',
  'input': 'msp-command-language-console-input',
  'keys': 'msp-command-language-console-keys',
  'kind': 'msp-command-language-console-kind',
  'label': 'msp-command-language-console-label',
  'line': 'msp-command-language-console-line',
  'lines': 'msp-command-language-console-lines',
  'mark': 'msp-command-language-console-mark',
  'notice': 'msp-command-language-console-notice',
  'noticeIcon': 'msp-command-language-console-noticeIcon',
  'output': 'msp-command-language-console-output',
  'picker': 'msp-command-language-console-picker',
  'segments': 'msp-command-language-console-segments',
  'send': 'msp-command-language-console-send',
  'side': 'msp-command-language-console-side',
  'sigActive': 'msp-command-language-console-sigActive',
  'sigName': 'msp-command-language-console-sigName',
  'signature': 'msp-command-language-console-signature',
  'spacer': 'msp-command-language-console-spacer',
  'spin': 'msp-command-language-console-spin',
  'status': 'msp-command-language-console-status',
  'summary': 'msp-command-language-console-summary',
  'swatch': 'msp-command-language-console-swatch',
  'time': 'msp-command-language-console-time',
  'typed': 'msp-command-language-console-typed',
  'usage': 'msp-command-language-console-usage',
  'value': 'msp-command-language-console-value',
} as const;

defineStyles('command-language/prism-viewer/console.module.css', `@layer msp.components {
.msp-command-language-console-console {
  flex-direction: column;
  flex: 1 1 0;
  height: 100%;
  min-height: 0;
  display: flex;
  position: relative;
}

.msp-command-language-console-bar {
  align-items: center;
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-2);
  display: flex;
}

.msp-command-language-console-segments {
  gap: var(--msp-space-1);
  display: inline-flex;
}

.msp-command-language-console-spacer {
  flex: 1;
}

.msp-command-language-console-picker {
  display: none;
}

.msp-command-language-console-output {
  gap: var(--msp-space-1);
  min-height: 0;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  scrollbar-gutter: stable;
  flex-direction: column;
  flex: 1;
  display: flex;
  overflow: auto;
}

.msp-command-language-console-entry {
  align-items: start;
  column-gap: var(--msp-space-2);
  row-gap: var(--msp-space-1);
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  grid-template-columns: 14px minmax(0, 1fr) auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-command-language-console-mark {
  margin-top: calc(var(--msp-space-1) / 2);
  color: var(--msp-color-success);
  display: inline-flex;

  .msp-command-language-console-entry[data-state="running"] & {
    color: var(--msp-color-text-subtle);
    animation: .9s linear infinite msp-command-language-console-spin;
  }

  .msp-command-language-console-entry[data-state="failed"] & {
    color: var(--msp-color-danger);
  }
}

@keyframes msp-command-language-console-spin {
  to {
    rotate: 360deg;
  }
}

.msp-command-language-console-command {
  align-items: baseline;
  gap: var(--msp-space-2);
  min-width: 0;
  display: flex;
}

.msp-command-language-console-chip {
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-well);
  color: var(--msp-color-text-subtle);
  letter-spacing: .02em;
  text-transform: uppercase;
  flex: none;
  font-size: 10px;
  font-weight: 600;
  line-height: 16px;
}

.msp-command-language-console-code {
  min-width: 0;
  color: var(--msp-color-text);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-sm);
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  line-height: 1.45;
}

.msp-command-language-console-side {
  align-items: center;
  gap: var(--msp-space-1);
  min-height: 20px;
  display: flex;
}

.msp-command-language-console-time {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
}

.msp-command-language-console-action {
  border-radius: var(--msp-radius-sm);
  width: 20px;
  height: 20px;
  color: var(--msp-color-text-muted);
  opacity: 0;
  cursor: pointer;
  transition: opacity var(--msp-motion-fast) var(--msp-ease);
  background: none;
  border: 0;
  justify-content: center;
  align-items: center;
  display: inline-flex;

  .msp-command-language-console-entry:hover &, &:focus-visible {
    opacity: 1;
  }

  @media (hover: none) {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-command-language-console-lines {
  gap: calc(var(--msp-space-1) / 2);
  grid-column: 2 / -1;
  min-width: 0;
  display: grid;
}

.msp-command-language-console-line {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  overflow-wrap: anywhere;
  line-height: 1.45;
}

.msp-command-language-console-value {
  max-width: 100%;
  padding: var(--msp-space-1) var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  color: var(--msp-color-text);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-sm);
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  justify-self: start;
}

.msp-command-language-console-notice {
  column-gap: var(--msp-space-1);
  font-size: var(--msp-font-size-sm);
  overflow-wrap: anywhere;
  grid-template-columns: 14px minmax(0, 1fr);
  align-items: start;
  line-height: 1.45;
  display: grid;

  &[data-kind="warning"] {
    color: var(--msp-color-warning);
  }

  &[data-kind="note"] {
    color: var(--msp-color-text-subtle);
  }

  &[data-kind="error"] {
    color: var(--msp-color-danger);
  }
}

.msp-command-language-console-noticeIcon {
  margin-top: calc(var(--msp-space-1) / 2);
  display: inline-flex;
}

.msp-command-language-console-hint {
  color: var(--msp-color-text-muted);
  grid-column: 2;
}

.msp-command-language-console-empty {
  gap: var(--msp-space-2);
  padding: var(--msp-space-3) var(--msp-space-2);
  display: grid;
}

.msp-command-language-console-emptyTitle {
  color: var(--msp-color-text);
  font-weight: 600;
}

.msp-command-language-console-emptyText {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  line-height: 1.5;
}

.msp-command-language-console-examples {
  gap: var(--msp-space-1);
  flex-wrap: wrap;
  display: flex;
}

.msp-command-language-console-example {
  max-width: 100%;
  padding: calc(var(--msp-space-1) / 2) var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  cursor: pointer;
  border: 0;
  line-height: 18px;
  overflow: hidden;

  &:hover {
    background: var(--msp-color-pressed);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 1px;
  }
}

.msp-command-language-console-composer {
  gap: var(--msp-space-1);
  margin: 0 var(--msp-space-2) var(--msp-space-2);
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-color-control);
  box-shadow: inset 0 0 0 1px var(--msp-field-edge);
  transition: box-shadow var(--msp-motion-fast) var(--msp-ease);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
  position: relative;

  &:focus-within {
    box-shadow: inset 0 0 0 1px var(--msp-field-edge),
      0 0 0 2px var(--msp-color-focus-ring);
  }
}

.msp-command-language-console-card {
  right: 0;
  bottom: calc(100% + var(--msp-space-3));
  z-index: var(--msp-z-overlay);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-color-panel);
  max-height: 280px;
  box-shadow: inset 0 0 0 1px var(--msp-field-edge),
    var(--msp-shadow-float);
  flex-direction: column;
  display: flex;
  position: absolute;
  left: 0;
}

.msp-command-language-console-signature {
  gap: var(--msp-space-1);
  padding: var(--msp-space-2);
  font-size: var(--msp-font-size-xs);
  flex: none;
  line-height: 1.4;
  display: grid;

  &:not(:last-child) {
    box-shadow: inset 0 -1px 0 var(--msp-field-edge);
  }
}

.msp-command-language-console-usage {
  color: var(--msp-color-text-muted);
  font-family: var(--msp-font-mono);
  overflow-wrap: anywhere;
}

.msp-command-language-console-summary {
  color: var(--msp-color-text-subtle);
}

.msp-command-language-console-completions {
  min-height: 0;
  padding: var(--msp-space-1);
  overflow-y: auto;
}

.msp-command-language-console-completion {
  align-items: center;
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  font-size: var(--msp-font-size-xs);
  cursor: pointer;
  grid-template-columns: 36px auto minmax(0, 1fr);
  line-height: 18px;
  display: grid;

  &[aria-selected="true"] {
    background: var(--msp-color-hover);
  }
}

.msp-command-language-console-kind {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  align-items: center;
  display: flex;
}

.msp-command-language-console-swatch {
  border-radius: var(--msp-radius-sm);
  width: 12px;
  height: 12px;
  box-shadow: inset 0 0 0 1px var(--msp-field-edge);
}

.msp-command-language-console-label {
  color: var(--msp-color-text);
  font-family: var(--msp-font-mono);
  white-space: nowrap;
}

.msp-command-language-console-typed {
  color: var(--msp-color-accent);
  background: none;
  font-weight: 600;
}

.msp-command-language-console-detail {
  color: var(--msp-color-text-muted);
  text-align: right;
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
}

.msp-command-language-console-field {
  align-items: flex-end;
  gap: var(--msp-space-2);
  display: flex;
}

.msp-command-language-console-input {
  min-width: 0;
  max-height: 160px;
  padding: calc(var(--msp-space-1) / 2) 0;
  color: var(--msp-color-text);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-sm);
  field-sizing: content;
  resize: none;
  background: none;
  border: 0;
  outline: 0;
  flex: 1;
  line-height: 1.45;

  &::placeholder {
    color: var(--msp-color-text-subtle);
  }
}

.msp-command-language-console-send {
  background: var(--msp-color-accent);
  width: 24px;
  height: 24px;
  color: var(--msp-color-on-accent);
  cursor: pointer;
  transition: opacity var(--msp-motion-fast) var(--msp-ease);
  border: 0;
  border-radius: 50%;
  flex: none;
  justify-content: center;
  align-items: center;
  display: inline-flex;

  &:disabled {
    background: var(--msp-color-well);
    color: var(--msp-color-text-subtle);
    cursor: default;
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 1px;
  }
}

.msp-command-language-console-status {
  align-items: center;
  gap: var(--msp-space-2);
  min-width: 0;
  height: 20px;
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-command-language-console-does {
  min-width: 0;
  color: var(--msp-color-text-muted);
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  overflow: hidden;

  &[data-kind="error"] {
    color: var(--msp-color-danger);
  }
}

.msp-command-language-console-hint {
  color: var(--msp-color-text-muted);
}

.msp-command-language-console-sigName {
  color: var(--msp-color-text);
}

.msp-command-language-console-sigActive {
  color: var(--msp-color-text);
  font-weight: 600;
}

.msp-command-language-console-keys {
  color: var(--msp-color-text-subtle);
  white-space: nowrap;
  flex: none;
}

.msp-command-language-console-drop {
  inset: var(--msp-space-2);
  padding: var(--msp-space-4);
  border: 1.5px dashed var(--msp-color-accent);
  border-radius: var(--msp-radius-lg);
  background: color-mix(in oklab, var(--msp-color-panel), transparent 8%);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  text-align: center;
  pointer-events: none;
  justify-content: center;
  align-items: center;
  display: none;
  position: absolute;

  .msp-command-language-console-console[data-dropping] & {
    display: flex;
  }
}
}
`);
