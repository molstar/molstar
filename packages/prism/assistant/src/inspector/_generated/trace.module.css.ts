/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from trace.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'at': 'msp-assistant-trace-at',
  'attempt': 'msp-assistant-trace-attempt',
  'attemptText': 'msp-assistant-trace-attemptText',
  'bar': 'msp-assistant-trace-bar',
  'bars': 'msp-assistant-trace-bars',
  'block': 'msp-assistant-trace-block',
  'column': 'msp-assistant-trace-column',
  'entries': 'msp-assistant-trace-entries',
  'entry': 'msp-assistant-trace-entry',
  'entryBody': 'msp-assistant-trace-entryBody',
  'entryRole': 'msp-assistant-trace-entryRole',
  'entryText': 'msp-assistant-trace-entryText',
  'entryTitle': 'msp-assistant-trace-entryTitle',
  'entryTokens': 'msp-assistant-trace-entryTokens',
  'fresh': 'msp-assistant-trace-fresh',
  'key': 'msp-assistant-trace-key',
  'lane': 'msp-assistant-trace-lane',
  'lanes': 'msp-assistant-trace-lanes',
  'legend': 'msp-assistant-trace-legend',
  'makeup': 'msp-assistant-trace-makeup',
  'milestones': 'msp-assistant-trace-milestones',
  'ordinal': 'msp-assistant-trace-ordinal',
  'overview': 'msp-assistant-trace-overview',
  'part': 'msp-assistant-trace-part',
  'pulse': 'msp-assistant-trace-pulse',
  'segment': 'msp-assistant-trace-segment',
  'swatch': 'msp-assistant-trace-swatch',
  'timing': 'msp-assistant-trace-timing',
  'track': 'msp-assistant-trace-track',
  'window': 'msp-assistant-trace-window',
} as const;

defineStyles('assistant/inspector/trace.module.css', `@layer msp.components {
.msp-assistant-trace-overview {
  align-items: center;
  gap: var(--msp-space-4);
  padding: var(--msp-space-4) var(--msp-space-5);
  flex: none;
  grid-template-columns: 40px minmax(0, 1fr) auto;
  display: grid;
}

.msp-assistant-trace-lanes {
  row-gap: var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  text-align: right;
  grid-template-rows: 12px 8px;
  font-size: 10px;
  line-height: 1;
  display: grid;

  & span:last-child {
    align-self: end;
  }
}

.msp-assistant-trace-track {
  gap: var(--msp-space-4);
  min-width: 0;
  display: flex;
}

.msp-assistant-trace-column {
  align-items: stretch;
  gap: var(--msp-space-1);
  min-width: 0;
  height: 22px;
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  flex-basis: 0;
  display: flex;

  &[data-selected] {
    background: var(--msp-color-control);
  }
}

.msp-assistant-trace-block {
  border-radius: var(--msp-radius-sm);
  min-width: 3px;
  transition: opacity var(--msp-motion-fast) var(--msp-ease);
  flex-basis: 0;
  padding: 0;

  &:hover {
    opacity: .7;
  }

  &[data-kind="request"] {
    background: var(--msp-color-text-subtle);
    flex: 0 0 3px;
    height: 100%;
  }

  &[data-kind="model"] {
    background: var(--msp-color-text-muted);
    align-self: flex-start;
    height: 12px;
  }

  &[data-kind="tool"] {
    background: var(--msp-color-border-strong);
    border-radius: 4px;
    align-self: flex-end;
    height: 8px;
  }

  &[data-nested] {
    opacity: .55;
  }

  &[data-failed] {
    background: var(--msp-color-danger);
  }

  &[data-retried] {
    background-image: repeating-linear-gradient(135deg, transparent 0 2px, var(--msp-color-bg) 2px 4px);
  }

  &[data-running] {
    animation: msp-assistant-trace-pulse 1.2s var(--msp-ease) infinite;
  }

  &[aria-pressed="true"] {
    outline: 1.5px solid var(--msp-color-accent);
    outline-offset: 1.5px;
  }
}

@keyframes msp-assistant-trace-pulse {
  50% {
    opacity: .35;
  }
}

.msp-assistant-trace-bars {
  gap: var(--msp-space-2);
  display: grid;
}

.msp-assistant-trace-bar {
  align-items: center;
  gap: var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  grid-template-columns: 16px minmax(0, 1fr) 96px;
  display: grid;

  & > :last-child {
    text-align: right;
  }
}

.msp-assistant-trace-milestones {
  gap: var(--msp-space-1);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  padding: 0;
  list-style: none;
  display: grid;
}

.msp-assistant-trace-at {
  min-width: 56px;
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-variant-numeric: tabular-nums;
  display: inline-block;
}

.msp-assistant-trace-window {
  gap: var(--msp-space-3);
  display: grid;
}

.msp-assistant-trace-makeup {
  border-radius: var(--msp-radius-sm);
  gap: 1px;
  height: 8px;
  display: flex;
  overflow: hidden;
}

.msp-assistant-trace-part {
  flex-basis: 0;
  min-width: 2px;
}

.msp-assistant-trace-swatch {
  width: 8px;
  height: 8px;
  margin-right: var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  display: inline-block;
}

.msp-assistant-trace-part, .msp-assistant-trace-swatch {
  &[data-part="guidance"] {
    background: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 78%);
  }

  &[data-part="tools"] {
    background: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 64%);
  }

  &[data-part="state"] {
    background: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 48%);
  }

  &[data-part="conversation"] {
    background: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 28%);
  }
}

.msp-assistant-trace-legend {
  gap: var(--msp-space-2) var(--msp-space-4);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  flex-wrap: wrap;
  display: flex;
}

.msp-assistant-trace-entries {
  gap: 1px;
  margin: 0;
  padding: 0;
  list-style: none;
  display: grid;
}

.msp-assistant-trace-entry {
  align-items: baseline;
  gap: var(--msp-space-3);
  width: 100%;
  padding: var(--msp-space-2) var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  text-align: left;
  grid-template-columns: 10px 64px minmax(0, auto) minmax(0, 1fr) auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-expanded="true"] {
    background: var(--msp-color-control);
  }
}

.msp-assistant-trace-fresh {
  border-radius: 50%;
  align-self: center;
  width: 5px;
  height: 5px;

  .msp-assistant-trace-entry[data-fresh] & {
    background: var(--msp-color-text-muted);
  }
}

.msp-assistant-trace-entryRole {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  letter-spacing: .04em;
  text-transform: uppercase;
  font-size: 10px;

  .msp-assistant-trace-entry[data-failed] & {
    color: var(--msp-color-danger);
  }
}

.msp-assistant-trace-entryTitle {
  white-space: nowrap;
  min-width: 0;
  overflow: hidden;
}

.msp-assistant-trace-entryText {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-assistant-trace-entryTokens {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
}

.msp-assistant-trace-entryBody {
  margin: var(--msp-space-1) 0 var(--msp-space-3);
}

.msp-assistant-trace-timing {
  gap: var(--msp-space-3);
  display: grid;
}

.msp-assistant-trace-attempt {
  align-items: center;
  gap: var(--msp-space-1) var(--msp-space-3);
  grid-template-columns: 24px minmax(0, 1fr);
  display: grid;
}

.msp-assistant-trace-ordinal {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
}

.msp-assistant-trace-lane {
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  height: 10px;
  position: relative;
}

.msp-assistant-trace-segment {
  border-radius: var(--msp-radius-sm);
  min-width: 2px;
  position: absolute;
  top: 0;
  bottom: 0;
}

.msp-assistant-trace-segment, .msp-assistant-trace-key i {
  &[data-kind="wait"] {
    background: repeating-linear-gradient(135deg, var(--msp-color-border-strong) 0 2px, transparent 2px 5px);
  }

  &[data-kind="pre"] {
    background: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 64%);
  }

  &[data-kind="out"] {
    background: var(--msp-color-text-muted);
  }

  &[data-kind="failed"] {
    background: var(--msp-color-danger);
  }

  &[data-kind="running"] {
    background: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 64%);
    animation: msp-assistant-trace-pulse 1.2s var(--msp-ease) infinite;
  }
}

.msp-assistant-trace-attemptText {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  grid-column: 2;

  &[data-failed] {
    color: var(--msp-color-danger);
  }
}

.msp-assistant-trace-key {
  gap: var(--msp-space-2) var(--msp-space-4);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex-wrap: wrap;
  display: flex;

  & i {
    width: 14px;
    height: 8px;
    margin-right: var(--msp-space-2);
    border-radius: var(--msp-radius-sm);
    display: inline-block;
  }
}

@container inspector (width <= 760px) {
  .msp-assistant-trace-overview {
    padding-inline: var(--msp-space-4);
    grid-template-columns: minmax(0, 1fr) auto;
  }

  .msp-assistant-trace-lanes {
    display: none;
  }
}

@media (prefers-reduced-motion: reduce) {
  .msp-assistant-trace-block, .msp-assistant-trace-segment {
    animation: none !important;
  }
}
}
`);
