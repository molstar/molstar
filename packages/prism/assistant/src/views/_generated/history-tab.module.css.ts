/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from history-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'answer': 'msp-assistant-history-tab-answer',
  'checkpoint': 'msp-assistant-history-tab-checkpoint',
  'checkpointMark': 'msp-assistant-history-tab-checkpointMark',
  'checkpointText': 'msp-assistant-history-tab-checkpointText',
  'chevron': 'msp-assistant-history-tab-chevron',
  'empty': 'msp-assistant-history-tab-empty',
  'fold': 'msp-assistant-history-tab-fold',
  'item': 'msp-assistant-history-tab-item',
  'list': 'msp-assistant-history-tab-list',
  'listContent': 'msp-assistant-history-tab-listContent',
  'prompt': 'msp-assistant-history-tab-prompt',
  'runs': 'msp-assistant-history-tab-runs',
  'spark': 'msp-assistant-history-tab-spark',
  'switch': 'msp-assistant-history-tab-switch',
  'switchButton': 'msp-assistant-history-tab-switchButton',
  'tab': 'msp-assistant-history-tab-tab',
  'turn': 'msp-assistant-history-tab-turn',
  'turnHead': 'msp-assistant-history-tab-turnHead',
  'undo': 'msp-assistant-history-tab-undo',
  'when': 'msp-assistant-history-tab-when',
  'who': 'msp-assistant-history-tab-who',
} as const;

defineStyles('assistant/views/history-tab.module.css', `@layer msp.components {
.msp-assistant-history-tab-tab {
  gap: var(--msp-space-2);
  height: 100%;
  min-height: 0;
  padding: 0 var(--msp-panel-inset) var(--msp-panel-inset);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  flex: 1 1 0;
  display: flex;
}

.msp-assistant-history-tab-switch {
  gap: var(--msp-space-1);
  flex: none;
  display: flex;
}

.msp-assistant-history-tab-switchButton {
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

.msp-assistant-history-tab-list {
  flex: 1 1 0;
}

.msp-assistant-history-tab-listContent {
  gap: var(--msp-space-1);
  display: grid;
}

.msp-assistant-history-tab-empty {
  padding: var(--msp-space-2) var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-assistant-history-tab-item {
  align-items: start;
  column-gap: var(--msp-space-2);
  padding: 0 var(--msp-inset);
  grid-template-columns: minmax(0, 1fr) auto;
  display: grid;

  & > .msp-assistant-history-tab-when {
    padding-top: 6px;
  }
}

.msp-assistant-history-tab-when {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
}

.msp-assistant-history-tab-turn {
  gap: var(--msp-space-1);
  margin: var(--msp-space-1) 0;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  box-shadow: inset 0 0 0 1px var(--msp-glass-edge);
  display: grid;

  &[data-undone] {
    opacity: .6;
  }
}

.msp-assistant-history-tab-turnHead {
  align-items: center;
  gap: var(--msp-space-2);
  display: flex;
}

.msp-assistant-history-tab-spark {
  color: var(--msp-color-accent);
  flex: none;
}

.msp-assistant-history-tab-who {
  min-width: 0;
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  flex: 1;
}

.msp-assistant-history-tab-prompt {
  margin: 0;
  font-weight: 600;
}

.msp-assistant-history-tab-answer {
  color: var(--msp-color-text-muted);
  margin: 0;
  line-height: 1.45;
}

.msp-assistant-history-tab-fold {
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

.msp-assistant-history-tab-chevron {
  transition: rotate var(--msp-motion-fast) var(--msp-ease);

  .msp-assistant-history-tab-fold[aria-expanded="true"] & {
    rotate: 180deg;
  }
}

.msp-assistant-history-tab-runs {
  display: grid;
}

.msp-assistant-history-tab-checkpoint {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 24px;
  margin-top: var(--msp-space-1);
  padding-top: var(--msp-space-1);
  box-shadow: inset 0 1px 0 var(--msp-glass-edge);
  display: flex;
}

.msp-assistant-history-tab-checkpointMark {
  color: var(--msp-color-text-subtle);
}

.msp-assistant-history-tab-checkpointText {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: 1;
}

.msp-assistant-history-tab-undo {
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
