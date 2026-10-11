/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from history-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'answer': 'msp-agent-history-tab-answer',
  'checkpoint': 'msp-agent-history-tab-checkpoint',
  'checkpointMark': 'msp-agent-history-tab-checkpointMark',
  'checkpointText': 'msp-agent-history-tab-checkpointText',
  'chevron': 'msp-agent-history-tab-chevron',
  'empty': 'msp-agent-history-tab-empty',
  'fold': 'msp-agent-history-tab-fold',
  'item': 'msp-agent-history-tab-item',
  'list': 'msp-agent-history-tab-list',
  'listContent': 'msp-agent-history-tab-listContent',
  'prompt': 'msp-agent-history-tab-prompt',
  'runs': 'msp-agent-history-tab-runs',
  'spark': 'msp-agent-history-tab-spark',
  'switch': 'msp-agent-history-tab-switch',
  'switchButton': 'msp-agent-history-tab-switchButton',
  'tab': 'msp-agent-history-tab-tab',
  'turn': 'msp-agent-history-tab-turn',
  'turnHead': 'msp-agent-history-tab-turnHead',
  'undo': 'msp-agent-history-tab-undo',
  'when': 'msp-agent-history-tab-when',
  'who': 'msp-agent-history-tab-who',
} as const;

defineStyles('views/agent/history-tab.module.css', `@layer msp.components {
.msp-agent-history-tab-tab {
  gap: var(--msp-space-2);
  height: 100%;
  min-height: 0;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  flex: 1 1 0;
  display: flex;
}

.msp-agent-history-tab-switch {
  gap: var(--msp-space-1);
  flex: none;
  display: flex;
}

.msp-agent-history-tab-switchButton {
  height: 24px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-selected="true"] {
    background: var(--msp-color-control);
    color: var(--msp-color-text);
  }
}

.msp-agent-history-tab-list {
  flex: 1 1 0;
}

.msp-agent-history-tab-listContent {
  gap: var(--msp-space-1);
  display: grid;
}

.msp-agent-history-tab-empty {
  padding: var(--msp-space-2) var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-agent-history-tab-item {
  align-items: start;
  column-gap: var(--msp-space-2);
  padding: 0 var(--msp-space-3);
  grid-template-columns: minmax(0, 1fr) auto;
  display: grid;

  & > .msp-agent-history-tab-when {
    padding-top: 6px;
  }
}

.msp-agent-history-tab-when {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
}

.msp-agent-history-tab-turn {
  gap: var(--msp-space-1);
  margin: var(--msp-space-1) 0;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  display: grid;

  &[data-undone] {
    opacity: .6;
  }
}

.msp-agent-history-tab-turnHead {
  align-items: center;
  gap: var(--msp-space-2);
  display: flex;
}

.msp-agent-history-tab-spark {
  color: var(--msp-color-accent);
  flex: none;
}

.msp-agent-history-tab-who {
  min-width: 0;
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  flex: 1;
}

.msp-agent-history-tab-prompt {
  margin: 0;
  font-weight: 600;
}

.msp-agent-history-tab-answer {
  color: var(--msp-color-text-muted);
  margin: 0;
  line-height: 1.45;
}

.msp-agent-history-tab-fold {
  align-items: center;
  gap: var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  justify-self: start;
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-agent-history-tab-chevron {
  transition: rotate var(--msp-motion-fast) var(--msp-ease);

  .msp-agent-history-tab-fold[aria-expanded="true"] & {
    rotate: 180deg;
  }
}

.msp-agent-history-tab-runs {
  display: grid;
}

.msp-agent-history-tab-checkpoint {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 24px;
  margin-top: var(--msp-space-1);
  padding-top: var(--msp-space-1);
  box-shadow: inset 0 1px 0 var(--msp-glass-edge);
  display: flex;
}

.msp-agent-history-tab-checkpointMark {
  color: var(--msp-color-text-subtle);
}

.msp-agent-history-tab-checkpointText {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: 1;
}

.msp-agent-history-tab-undo {
  height: 22px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-accent-strong);
  font-size: var(--msp-font-size-xs);
  font-weight: 600;

  &:hover {
    background: var(--msp-color-hover);
  }
}
}
`);
