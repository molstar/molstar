/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from agent-turn.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'allow': 'msp-agent-agent-turn-allow',
  'approval': 'msp-agent-agent-turn-approval',
  'approvalText': 'msp-agent-agent-turn-approvalText',
  'approvalWhy': 'msp-agent-agent-turn-approvalWhy',
  'bar': 'msp-agent-agent-turn-bar',
  'barFill': 'msp-agent-agent-turn-barFill',
  'body': 'msp-agent-agent-turn-body',
  'bodyContent': 'msp-agent-agent-turn-bodyContent',
  'card': 'msp-agent-agent-turn-card',
  'deny': 'msp-agent-agent-turn-deny',
  'footer': 'msp-agent-agent-turn-footer',
  'header': 'msp-agent-agent-turn-header',
  'icon': 'msp-agent-agent-turn-icon',
  'line': 'msp-agent-agent-turn-line',
  'mark': 'msp-agent-agent-turn-mark',
  'note': 'msp-agent-agent-turn-note',
  'prompt': 'msp-agent-agent-turn-prompt',
  'run': 'msp-agent-agent-turn-run',
  'runDetail': 'msp-agent-agent-turn-runDetail',
  'runTime': 'msp-agent-agent-turn-runTime',
  'runTitle': 'msp-agent-agent-turn-runTitle',
  'runs': 'msp-agent-agent-turn-runs',
  'spark': 'msp-agent-agent-turn-spark',
  'spark-pulse': 'msp-agent-agent-turn-spark-pulse',
  'spin': 'msp-agent-agent-turn-spin',
  'state': 'msp-agent-agent-turn-state',
  'stop': 'msp-agent-agent-turn-stop',
  'summary': 'msp-agent-agent-turn-summary',
  'text': 'msp-agent-agent-turn-text',
  'undo': 'msp-agent-agent-turn-undo',
  'undone': 'msp-agent-agent-turn-undone',
} as const;

defineStyles('agent/views/agent-turn.module.css', `@layer msp.components {
.msp-agent-agent-turn-card {
  margin-top: var(--msp-space-1);
  width: min(380px, calc(100cqw - 2 * var(--msp-space-4)));
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  max-height: min(420px, 60cqh);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  display: flex;
}

.msp-agent-agent-turn-header {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 34px;
  padding: 0 var(--msp-space-2) 0 var(--msp-space-3);
  display: flex;
}

.msp-agent-agent-turn-spark {
  color: var(--msp-color-accent);
  flex: none;

  .msp-agent-agent-turn-card:is([data-status="thinking"], [data-status="working"]) & {
    animation: msp-agent-agent-turn-spark-pulse 1.4s var(--msp-ease) infinite;
  }
}

@keyframes msp-agent-agent-turn-spark-pulse {
  50% {
    opacity: .35;
  }
}

.msp-agent-agent-turn-prompt {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  font-weight: 600;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-agent-agent-turn-state {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;

  .msp-agent-agent-turn-card[data-status="waiting"] & {
    color: var(--msp-color-accent);
  }

  .msp-agent-agent-turn-card[data-status="failed"] & {
    color: var(--msp-color-danger);
  }
}

.msp-agent-agent-turn-icon, .msp-agent-agent-turn-stop {
  justify-content: center;
  align-items: center;
  gap: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  height: 24px;
  color: var(--msp-color-text-muted);
  flex: none;
  display: inline-flex;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-agent-agent-turn-icon {
  width: 24px;
}

.msp-agent-agent-turn-stop {
  padding: 0 var(--msp-space-2);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-xs);
}

.msp-agent-agent-turn-body {
  flex: 0 auto;
}

.msp-agent-agent-turn-bodyContent {
  gap: var(--msp-space-2);
  padding: 0 var(--msp-space-3) var(--msp-space-2);
  display: grid;
}

.msp-agent-agent-turn-text {
  margin: 0;
  line-height: 1.45;

  &[data-streaming] {
    color: var(--msp-color-text-muted);
  }
}

.msp-agent-agent-turn-note {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-agent-agent-turn-runs {
  display: grid;
}

.msp-agent-agent-turn-run {
  align-items: center;
  column-gap: var(--msp-space-2);
  min-height: 24px;
  padding: var(--msp-space-1) 0;
  grid-template-columns: 14px minmax(0, 1fr) auto;
  display: grid;

  &[data-read] {
    color: var(--msp-color-text-muted);
  }
}

.msp-agent-agent-turn-mark {
  color: var(--msp-color-text-subtle);

  .msp-agent-agent-turn-run[data-status="running"] & {
    color: var(--msp-color-text);
    animation: .9s linear infinite msp-agent-agent-turn-spin;
  }

  .msp-agent-agent-turn-run[data-status="waiting"] & {
    color: var(--msp-color-accent);
  }

  .msp-agent-agent-turn-run[data-status="failed"] & {
    color: var(--msp-color-danger);
  }
}

@keyframes msp-agent-agent-turn-spin {
  to {
    rotate: 360deg;
  }
}

.msp-agent-agent-turn-runTitle {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;

  .msp-agent-agent-turn-run[data-status="stopped"] &, .msp-agent-agent-turn-run[data-status="denied"] & {
    color: var(--msp-color-text-subtle);
    text-decoration: line-through;
  }
}

.msp-agent-agent-turn-runTime {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
}

.msp-agent-agent-turn-runDetail {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  grid-column: 2 / -1;

  .msp-agent-agent-turn-run[data-status="failed"] & {
    color: var(--msp-color-danger);
  }
}

.msp-agent-agent-turn-bar {
  height: 2px;
  margin-top: var(--msp-space-1);
  background: var(--msp-color-control);
  border-radius: 1px;
  grid-column: 2 / -1;
  overflow: hidden;
}

.msp-agent-agent-turn-barFill {
  background: var(--msp-color-text-muted);
  height: 100%;
  transition: width var(--msp-motion-fast) var(--msp-ease);
  display: block;
}

.msp-agent-agent-turn-approval {
  align-items: center;
  gap: var(--msp-space-2);
  margin: 0 var(--msp-space-2) var(--msp-space-2);
  padding: var(--msp-space-2) var(--msp-space-2) var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-accent), transparent 88%);
  display: flex;
}

.msp-agent-agent-turn-approvalText {
  flex: 1;
  min-width: 0;
}

.msp-agent-agent-turn-approvalWhy {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  display: block;
}

.msp-agent-agent-turn-deny, .msp-agent-agent-turn-allow {
  height: 24px;
  padding: 0 var(--msp-inset);
  border-radius: var(--msp-radius-md);
  flex: none;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-agent-agent-turn-deny {
  color: var(--msp-color-text-muted);
}

.msp-agent-agent-turn-allow {
  background: color-mix(in oklab, var(--msp-color-accent), transparent 80%);
  color: var(--msp-color-accent-strong);
  font-weight: 600;

  &:hover {
    background: color-mix(in oklab, var(--msp-color-accent), transparent 70%);
  }
}

.msp-agent-agent-turn-footer {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 32px;
  padding: 0 var(--msp-space-2) 0 var(--msp-space-3);
  box-shadow: inset 0 1px 0 var(--msp-glass-edge);
  display: flex;
}

.msp-agent-agent-turn-summary {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-agent-agent-turn-undone {
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
}

.msp-agent-agent-turn-undo {
  align-items: center;
  gap: var(--msp-space-1);
  height: 24px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-accent-strong);
  font-weight: 600;
  display: inline-flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-agent-agent-turn-line {
  display: contents;

  &:is([data-status="thinking"], [data-status="working"]) .msp-agent-agent-turn-spark {
    animation: msp-agent-agent-turn-spark-pulse 1.4s var(--msp-ease) infinite;
  }

  &[data-status="waiting"] .msp-agent-agent-turn-state {
    color: var(--msp-color-accent);
  }
}
}
`);
