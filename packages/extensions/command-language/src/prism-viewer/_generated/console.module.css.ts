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
  'console': 'msp-command-language-console-console',
  'dialects': 'msp-command-language-console-dialects',
  'entry': 'msp-command-language-console-entry',
  'error': 'msp-command-language-console-error',
  'field': 'msp-command-language-console-field',
  'foot': 'msp-command-language-console-foot',
  'hint': 'msp-command-language-console-hint',
  'input': 'msp-command-language-console-input',
  'line': 'msp-command-language-console-line',
  'mono': 'msp-command-language-console-mono',
  'output': 'msp-command-language-console-output',
  'picker': 'msp-command-language-console-picker',
  'prompt': 'msp-command-language-console-prompt',
  'quiet': 'msp-command-language-console-quiet',
  'spacer': 'msp-command-language-console-spacer',
  'toolbar': 'msp-command-language-console-toolbar',
  'warning': 'msp-command-language-console-warning',
} as const;

defineStyles('command-language/prism-viewer/console.module.css', `@layer msp.components {
.msp-command-language-console-console {
  gap: var(--msp-space-2);
  flex-direction: column;
  flex: 1;
  min-height: 0;
  display: flex;
}

.msp-command-language-console-toolbar {
  align-items: center;
  gap: var(--msp-space-1);
  padding: var(--msp-space-1) var(--msp-space-2) 0;
  display: flex;
}

.msp-command-language-console-dialects {
  gap: var(--msp-space-1);
  display: flex;
}

.msp-command-language-console-spacer {
  flex: 1;
}

.msp-command-language-console-picker {
  display: none;
}

.msp-command-language-console-output {
  min-height: 0;
  padding: 0 var(--msp-space-2);
  flex: 1;
  overflow: auto;
}

.msp-command-language-console-entry {
  gap: var(--msp-space-1);
  padding: var(--msp-space-1) 0;
  display: grid;
}

.msp-command-language-console-mono {
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-sm);
  white-space: pre-wrap;
  overflow-wrap: anywhere;
}

.msp-command-language-console-input {
  gap: var(--msp-space-2);
  display: flex;
}

.msp-command-language-console-prompt {
  color: var(--msp-color-text-muted);
  white-space: nowrap;
  flex: none;
}

.msp-command-language-console-line {
  color: var(--msp-color-text);
}

.msp-command-language-console-quiet {
  color: var(--msp-color-text-muted);
}

.msp-command-language-console-warning {
  color: var(--msp-color-warning);
}

.msp-command-language-console-error {
  color: var(--msp-color-danger);
}

.msp-command-language-console-hint {
  padding: var(--msp-space-2) 0;
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
}

.msp-command-language-console-foot {
  padding: 0 var(--msp-space-2) var(--msp-space-2);
}

.msp-command-language-console-field {
  box-sizing: border-box;
  width: 100%;
  min-height: 56px;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  box-shadow: inset 0 0 0 1px var(--msp-field-edge);
  color: var(--msp-color-text);
  resize: vertical;
  border: 0;
}
}
`);
