/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from parts.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'budget': 'msp-assistant-parts-budget',
  'budgetBar': 'msp-assistant-parts-budgetBar',
  'budgetText': 'msp-assistant-parts-budgetText',
  'code': 'msp-assistant-parts-code',
  'count': 'msp-assistant-parts-count',
  'dim': 'msp-assistant-parts-dim',
  'fact': 'msp-assistant-parts-fact',
  'facts': 'msp-assistant-parts-facts',
  'h': 'msp-assistant-parts-h',
  'head': 'msp-assistant-parts-head',
  'jKey': 'msp-assistant-parts-jKey',
  'jLit': 'msp-assistant-parts-jLit',
  'jStr': 'msp-assistant-parts-jStr',
  'json': 'msp-assistant-parts-json',
  'kind': 'msp-assistant-parts-kind',
  'link': 'msp-assistant-parts-link',
  'meter': 'msp-assistant-parts-meter',
  'meterFill': 'msp-assistant-parts-meterFill',
  'mono': 'msp-assistant-parts-mono',
  'note': 'msp-assistant-parts-note',
  'prose': 'msp-assistant-parts-prose',
  'tab': 'msp-assistant-parts-tab',
  'tabs': 'msp-assistant-parts-tabs',
  'usage': 'msp-assistant-parts-usage',
  'where': 'msp-assistant-parts-where',
} as const;

defineStyles('assistant/inspector/parts.module.css', `@layer msp.components {
.msp-assistant-parts-tabs {
  align-items: center;
  gap: var(--msp-space-1);
  flex: none;
  display: flex;
}

.msp-assistant-parts-tab {
  align-items: baseline;
  gap: var(--msp-space-2);
  height: 24px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  flex: none;
  line-height: 24px;
  display: inline-flex;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-selected="true"] {
    background: var(--msp-color-control);
    color: var(--msp-color-text);
  }
}

.msp-assistant-parts-count {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;

  &:empty {
    display: none;
  }
}

.msp-assistant-parts-head {
  align-items: baseline;
  gap: var(--msp-space-3);
  min-width: 0;
  display: flex;
}

.msp-assistant-parts-kind {
  color: var(--msp-color-text-muted);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  letter-spacing: .04em;
  text-transform: uppercase;
  flex: none;
}

.msp-assistant-parts-where {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-assistant-parts-facts {
  gap: var(--msp-space-1);
  margin: 0;
  display: grid;
}

.msp-assistant-parts-fact {
  gap: var(--msp-space-3);
  grid-template-columns: 96px minmax(0, 1fr);
  line-height: 1.45;
  display: grid;

  & dt {
    color: var(--msp-color-text-subtle);
  }

  & dd {
    overflow-wrap: anywhere;
    min-width: 0;
    margin: 0;
  }
}

.msp-assistant-parts-code {
  padding: var(--msp-space-3) var(--msp-space-4);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  color: var(--msp-color-text);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  white-space: pre-wrap;
  overflow-wrap: anywhere;
  margin: 0;
  line-height: 1.55;

  &[data-failed] {
    color: var(--msp-color-danger);
  }
}

.msp-assistant-parts-json {
  tab-size: 2;
}

.msp-assistant-parts-jKey {
  color: var(--msp-color-text-muted);
}

.msp-assistant-parts-jStr {
  color: color-mix(in oklab, var(--msp-color-success), var(--msp-color-text) 35%);
}

.msp-assistant-parts-jLit {
  color: color-mix(in oklab, var(--msp-color-info), var(--msp-color-text) 25%);
}

.msp-assistant-parts-usage {
  color: var(--msp-color-text-subtle);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}

.msp-assistant-parts-meter {
  background: var(--msp-color-control);
  border-radius: 2px;
  height: 3px;
  display: block;
  overflow: hidden;
}

.msp-assistant-parts-meterFill {
  border-radius: inherit;
  background: var(--msp-color-text-subtle);
  height: 100%;
  display: block;

  .msp-assistant-parts-meter[data-danger] & {
    background: var(--msp-color-danger);
  }
}

.msp-assistant-parts-budget {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  display: inline-flex;
}

.msp-assistant-parts-budgetBar {
  flex: none;
  width: 64px;

  @container (width <= 640px) {
    display: none;
  }
}

.msp-assistant-parts-budgetText {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-assistant-parts-note {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  line-height: 1.5;
}

.msp-assistant-parts-h {
  margin: var(--msp-space-2) 0 0;
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  font-weight: 600;
}

.msp-assistant-parts-prose {
  font-size: var(--msp-font-size-md);
  overflow-wrap: anywhere;
  margin: 0;
  line-height: 1.45;
}

.msp-assistant-parts-link {
  color: var(--msp-color-accent-strong);
  text-align: left;
  text-underline-offset: 2px;
  padding: 0;
  text-decoration: underline #0000;

  &:hover {
    text-decoration-color: currentColor;
  }
}

.msp-assistant-parts-mono {
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
}

.msp-assistant-parts-dim {
  color: var(--msp-color-text-subtle);
}
}
`);
