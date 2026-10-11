/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from assistant-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'content': 'msp-agent-assistant-tab-content',
  'hint': 'msp-agent-assistant-tab-hint',
  'tab': 'msp-agent-assistant-tab-tab',
} as const;

defineStyles('views/agent/assistant-tab.module.css', `@layer msp.components {
.msp-agent-assistant-tab-tab {
  flex: 1 1 0;
  height: 100%;
  min-height: 0;
}

.msp-agent-assistant-tab-content {
  align-content: start;
  gap: var(--msp-space-4);
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  display: grid;
}

.msp-agent-assistant-tab-hint {
  padding: 0 var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  line-height: 1.45;
}
}
`);
