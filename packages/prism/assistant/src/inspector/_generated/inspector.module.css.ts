/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from inspector.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'args': 'msp-assistant-inspector-args',
  'body': 'msp-assistant-inspector-body',
  'chip': 'msp-assistant-inspector-chip',
  'chips': 'msp-assistant-inspector-chips',
  'detail': 'msp-assistant-inspector-detail',
  'detailContent': 'msp-assistant-inspector-detailContent',
  'dot': 'msp-assistant-inspector-dot',
  'empty': 'msp-assistant-inspector-empty',
  'flag': 'msp-assistant-inspector-flag',
  'focus': 'msp-assistant-inspector-focus',
  'head': 'msp-assistant-inspector-head',
  'headMeta': 'msp-assistant-inspector-headMeta',
  'headTitle': 'msp-assistant-inspector-headTitle',
  'inspector': 'msp-assistant-inspector-inspector',
  'list': 'msp-assistant-inspector-list',
  'listContent': 'msp-assistant-inspector-listContent',
  'meta': 'msp-assistant-inspector-meta',
  'name': 'msp-assistant-inspector-name',
  'nest': 'msp-assistant-inspector-nest',
  'pulse': 'msp-assistant-inspector-pulse',
  'result': 'msp-assistant-inspector-result',
  'role': 'msp-assistant-inspector-role',
  'row': 'msp-assistant-inspector-row',
  'spacer': 'msp-assistant-inspector-spacer',
  'split': 'msp-assistant-inspector-split',
  'state': 'msp-assistant-inspector-state',
  'status': 'msp-assistant-inspector-status',
  'time': 'msp-assistant-inspector-time',
  'tokens': 'msp-assistant-inspector-tokens',
  'view': 'msp-assistant-inspector-view',
} as const;

defineStyles('assistant/inspector/inspector.module.css', `@layer msp.components {
.msp-assistant-inspector-inspector {
  min-height: 0;
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  flex: 1 1 0;
  display: flex;
  container: inspector / inline-size;
}

.msp-assistant-inspector-view {
  flex-direction: column;
  flex: 1 1 0;
  min-height: 0;
  display: flex;
}

.msp-assistant-inspector-split {
  flex: 1 1 0;
  grid-template-columns: minmax(0, 1fr) 380px;
  min-height: 0;
  display: grid;
}

.msp-assistant-inspector-list {
  height: 100%;
  min-height: 0;
}

.msp-assistant-inspector-listContent {
  padding-bottom: var(--msp-space-6);
}

.msp-assistant-inspector-detail {
  background: var(--msp-color-panel);
  height: 100%;
  min-height: 0;
}

.msp-assistant-inspector-detailContent {
  align-content: start;
  gap: var(--msp-space-4);
  min-width: 0;
  padding: var(--msp-space-5) var(--msp-space-5) var(--msp-space-6);
  display: grid;
}

.msp-assistant-inspector-head {
  z-index: 1;
  align-items: baseline;
  gap: var(--msp-space-3);
  width: 100%;
  min-width: 0;
  padding: var(--msp-space-4) var(--msp-space-5) var(--msp-space-2);
  background: var(--msp-color-bg);
  text-align: left;
  display: flex;
  position: sticky;
  top: 0;

  &:hover .msp-assistant-inspector-headTitle {
    text-decoration: underline;
    text-decoration-color: var(--msp-color-border-strong);
    text-underline-offset: 3px;
  }

  &[aria-selected="true"] {
    background: var(--msp-color-control);
  }
}

.msp-assistant-inspector-headTitle {
  flex: none;
  font-weight: 600;
}

.msp-assistant-inspector-headMeta {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  text-align: right;
  white-space: nowrap;
  flex: auto;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-assistant-inspector-status {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;

  &[data-status="failed"] {
    color: var(--msp-color-danger);
  }

  &[data-status="waiting"] {
    color: var(--msp-color-accent);
  }

  &:is([data-status="thinking"], [data-status="working"]) {
    color: var(--msp-color-text-muted);
    animation: msp-assistant-inspector-pulse 1.4s var(--msp-ease) infinite;
  }
}

@keyframes msp-assistant-inspector-pulse {
  50% {
    opacity: .4;
  }
}

.msp-assistant-inspector-row {
  align-items: baseline;
  gap: var(--msp-space-3);
  width: 100%;
  min-height: 28px;
  padding: var(--msp-space-2) var(--msp-space-5);
  text-align: left;
  grid-template-columns: 56px minmax(0, 1fr) auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-selected="true"] {
    background: var(--msp-color-control);
  }

  &[data-off] {
    color: var(--msp-color-text-subtle);
  }
}

.msp-assistant-inspector-role {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  letter-spacing: .04em;
  text-align: right;
  text-transform: uppercase;
  white-space: nowrap;
  font-size: 10px;
  overflow: hidden;

  .msp-assistant-inspector-row[data-failed] & {
    color: var(--msp-color-danger);
  }
}

.msp-assistant-inspector-body {
  align-items: baseline;
  gap: var(--msp-space-3);
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  display: flex;
  overflow: hidden;

  & > * {
    flex: none;
  }

  .msp-assistant-inspector-row[data-nested] & {
    padding-left: var(--msp-space-5);
  }
}

.msp-assistant-inspector-nest {
  color: var(--msp-color-text-subtle);
}

.msp-assistant-inspector-meta {
  align-items: baseline;
  gap: var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
  display: inline-flex;
}

.msp-assistant-inspector-time {
  text-align: right;
  min-width: 52px;
}

.msp-assistant-inspector-tokens {
  text-align: right;
  min-width: 80px;
}

.msp-assistant-inspector-flag {
  color: var(--msp-color-text-muted);
}

.msp-assistant-inspector-name {
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-weight: 600;
}

.msp-assistant-inspector-args {
  color: var(--msp-color-text-muted);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
}

.msp-assistant-inspector-result {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);

  .msp-assistant-inspector-row[data-failed] & {
    color: var(--msp-color-danger);
  }
}

.msp-assistant-inspector-focus {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-assistant-inspector-chips {
  gap: var(--msp-space-2);
  flex-wrap: wrap;
  display: inline-flex;
}

.msp-assistant-inspector-chip {
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  color: var(--msp-color-text-muted);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
}

.msp-assistant-inspector-empty {
  padding: var(--msp-space-6) var(--msp-space-5);
  color: var(--msp-color-text-subtle);
  margin: 0;
}

.msp-assistant-inspector-spacer {
  flex: 1 1 0;
  min-width: 0;
}

.msp-assistant-inspector-state {
  align-items: center;
  gap: var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
  display: inline-flex;
}

.msp-assistant-inspector-dot {
  background: var(--msp-color-border-strong);
  border-radius: 50%;
  width: 6px;
  height: 6px;

  .msp-assistant-inspector-state[data-state="working"] & {
    background: var(--msp-color-text-muted);
    animation: msp-assistant-inspector-pulse 1.4s var(--msp-ease) infinite;
  }

  .msp-assistant-inspector-state[data-state="waiting"] & {
    background: var(--msp-color-accent);
  }
}

@container (width <= 640px) {
  .msp-assistant-inspector-state {
    display: none;
  }
}

@media (prefers-reduced-motion: reduce) {
  .msp-assistant-inspector-status, .msp-assistant-inspector-dot {
    animation: none !important;
  }
}

@container inspector (width <= 760px) {
  .msp-assistant-inspector-split {
    grid-template-rows: minmax(0, 1fr) minmax(0, 46%);
    grid-template-columns: minmax(0, 1fr);
  }

  .msp-assistant-inspector-row {
    padding-inline: var(--msp-space-4);
    grid-template-columns: 44px minmax(0, 1fr) auto;
  }

  .msp-assistant-inspector-head {
    padding-inline: var(--msp-space-4);
  }

  .msp-assistant-inspector-detailContent {
    padding: var(--msp-space-4);
  }
}
}
`);
