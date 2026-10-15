/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from history.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'bar': 'msp-viewer-history-bar',
  'caption': 'msp-viewer-history-caption',
  'count': 'msp-viewer-history-count',
  'detail': 'msp-viewer-history-detail',
  'empty': 'msp-viewer-history-empty',
  'half': 'msp-viewer-history-half',
  'halfBody': 'msp-viewer-history-halfBody',
  'halfTitle': 'msp-viewer-history-halfTitle',
  'history': 'msp-viewer-history-history',
  'keep': 'msp-viewer-history-keep',
  'label': 'msp-viewer-history-label',
  'list': 'msp-viewer-history-list',
  'more': 'msp-viewer-history-more',
  'origin': 'msp-viewer-history-origin',
  'prompt': 'msp-viewer-history-prompt',
  'reply': 'msp-viewer-history-reply',
  'replyFull': 'msp-viewer-history-replyFull',
  'runs': 'msp-viewer-history-runs',
  'separator': 'msp-viewer-history-separator',
  'separatorText': 'msp-viewer-history-separatorText',
  'split': 'msp-viewer-history-split',
  'step': 'msp-viewer-history-step',
  'tabs': 'msp-viewer-history-tabs',
  'time': 'msp-viewer-history-time',
  'turn': 'msp-viewer-history-turn',
  'turnBody': 'msp-viewer-history-turnBody',
  'turnHead': 'msp-viewer-history-turnHead',
  'unsafe': 'msp-viewer-history-unsafe',
} as const;

defineStyles('viewer/ui/history.module.css', `@layer msp.components {
.msp-viewer-history-history {
  contain: content;
  gap: var(--msp-space-2);
  padding: 0 var(--msp-panel-inset) var(--msp-panel-inset);
  font-size: var(--msp-font-size-sm);
  display: grid;
}

.msp-viewer-history-list {
  contain: content;
  gap: var(--msp-space-2);
  display: grid;
}

.msp-viewer-history-bar {
  align-items: center;
  gap: var(--msp-space-3);
  display: flex;
}

.msp-viewer-history-count {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin-left: auto;
}

.msp-viewer-history-empty {
  padding: var(--msp-space-4);
  color: var(--msp-color-text-muted);
}

.msp-viewer-history-step, .msp-viewer-history-turn, .msp-viewer-history-separator {
  content-visibility: auto;
  contain-intrinsic-size: auto 2.25rem;
}

.msp-viewer-history-step {
  align-items: baseline;
  gap: var(--msp-space-1) var(--msp-space-3);
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  grid-template-columns: auto minmax(0, 1fr) auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[data-unsafe] .msp-viewer-history-label {
    color: var(--msp-color-text-muted);
  }
}

.msp-viewer-history-origin {
  color: var(--msp-color-text-subtle);
  align-self: center;
}

.msp-viewer-history-label {
  overflow-wrap: anywhere;
}

.msp-viewer-history-time {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
}

.msp-viewer-history-detail {
  align-items: center;
  gap: var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  grid-column: 2 / -1;
  display: flex;
}

.msp-viewer-history-unsafe {
  overflow-wrap: anywhere;
  opacity: .8;
  min-width: 0;
}

.msp-viewer-history-keep {
  min-height: var(--msp-space-5);
  align-items: center;
  display: flex;
}

.msp-viewer-history-caption {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-history-more {
  justify-content: center;
  display: flex;
}

.msp-viewer-history-separator {
  align-items: baseline;
  gap: var(--msp-space-3);
  padding: var(--msp-space-1) var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-viewer-history-separatorText {
  overflow-wrap: anywhere;
  flex: 1;
  min-width: 0;
}

.msp-viewer-history-turn {
  gap: var(--msp-space-1);
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-history-turnHead {
  align-items: baseline;
  gap: var(--msp-space-3);
  color: inherit;
  font: inherit;
  text-align: left;
  cursor: pointer;
  background: none;
  border: 0;
  grid-template-columns: auto minmax(0, 1fr) auto;
  padding: 0;
  display: grid;
}

.msp-viewer-history-reply, .msp-viewer-history-replyFull, .msp-viewer-history-prompt {
  color: var(--msp-color-text-muted);
  overflow-wrap: anywhere;
  margin: 0;
}

.msp-viewer-history-reply {
  max-height: 1.5em;
  overflow: hidden;
}

.msp-viewer-history-replyFull, .msp-viewer-history-prompt {
  white-space: pre-wrap;
}

.msp-viewer-history-turnBody {
  gap: var(--msp-space-2);
  display: grid;
}

.msp-viewer-history-runs {
  gap: var(--msp-space-1);
  padding-left: var(--msp-space-3);
  display: grid;
}

.msp-viewer-history-step[data-undone] .msp-viewer-history-label, .msp-viewer-history-turn[data-undone] .msp-viewer-history-label {
  color: var(--msp-color-text-muted);
  text-decoration: line-through;
}

.msp-viewer-history-step[data-failed] .msp-viewer-history-origin {
  color: var(--msp-color-danger);
}

.msp-viewer-history-step[data-nested] {
  padding-inline: 0;
}

.msp-viewer-history-split {
  gap: var(--msp-space-3);
  flex: 1 1 0;
  grid-template-rows: minmax(0, 1fr) minmax(0, 1fr);
  height: 100%;
  min-height: 0;
  display: grid;
}

.msp-viewer-history-half {
  flex-direction: column;
  min-height: 0;
  display: flex;
}

.msp-viewer-history-halfTitle {
  padding: var(--msp-space-2) calc(var(--msp-panel-inset) + var(--msp-space-3));
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  text-transform: uppercase;
  letter-spacing: .04em;
  font-weight: 600;
}

.msp-viewer-history-halfBody {
  flex: 1 1 0;
  min-height: 0;
}

.msp-viewer-history-tabs {
  gap: var(--msp-space-2);
  padding-inline: var(--msp-panel-inset);
  display: grid;
}
}
`);
