/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sequence-inspector.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'box': 'msp-viewer-sequence-inspector-box',
  'boxButton': 'msp-viewer-sequence-inspector-boxButton',
  'branched': 'msp-viewer-sequence-inspector-branched',
  'chevron': 'msp-viewer-sequence-inspector-chevron',
  'chevronOpen': 'msp-viewer-sequence-inspector-chevronOpen',
  'chip': 'msp-viewer-sequence-inspector-chip',
  'chips': 'msp-viewer-sequence-inspector-chips',
  'dot': 'msp-viewer-sequence-inspector-dot',
  'entities': 'msp-viewer-sequence-inspector-entities',
  'entity': 'msp-viewer-sequence-inspector-entity',
  'ligand': 'msp-viewer-sequence-inspector-ligand',
  'line': 'msp-viewer-sequence-inspector-line',
  'lineText': 'msp-viewer-sequence-inspector-lineText',
  'list': 'msp-viewer-sequence-inspector-list',
  'meta': 'msp-viewer-sequence-inspector-meta',
  'more': 'msp-viewer-sequence-inspector-more',
  'moreCount': 'msp-viewer-sequence-inspector-moreCount',
  'name': 'msp-viewer-sequence-inspector-name',
  'polymer': 'msp-viewer-sequence-inspector-polymer',
  'row': 'msp-viewer-sequence-inspector-row',
  'status': 'msp-viewer-sequence-inspector-status',
  'statusText': 'msp-viewer-sequence-inspector-statusText',
  'text': 'msp-viewer-sequence-inspector-text',
  'word': 'msp-viewer-sequence-inspector-word',
} as const;

defineStyles('viewer/views/structure/sequence-inspector.module.css', `@layer msp.components {
.msp-viewer-sequence-inspector-status {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 28px;
  display: flex;
}

.msp-viewer-sequence-inspector-statusText {
  min-width: 0;
  color: var(--msp-color-text-muted);
  white-space: nowrap;
  flex: 1;
  overflow: hidden;
}

.msp-viewer-sequence-inspector-word {
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-sequence-inspector-line {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 28px;
  color: var(--msp-color-text-muted);
  text-align: left;
  display: flex;
}

.msp-viewer-sequence-inspector-lineText {
  min-width: 0;
}

.msp-viewer-sequence-inspector-entities {
  display: grid;
}

.msp-viewer-sequence-inspector-list {
  margin: 0 calc(-1 * var(--msp-space-2));
  display: grid;
}

.msp-viewer-sequence-inspector-entity {
  min-width: 0;
  display: grid;
}

.msp-viewer-sequence-inspector-row {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 36px;
  padding: 0 var(--msp-space-2) 0 0;
  border-radius: var(--msp-radius-md);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-sequence-inspector-box {
  border: 1px solid var(--msp-color-border-strong);
  border-radius: var(--msp-radius-sm);
  width: 16px;
  height: 16px;
  color: var(--msp-color-on-accent);
  flex: none;
  place-items: center;
  display: grid;

  &[data-checked="true"], &[data-checked="mixed"] {
    border-color: var(--msp-color-accent);
    background: var(--msp-color-accent);
  }
}

.msp-viewer-sequence-inspector-boxButton {
  flex: none;
  place-items: center;
  width: 24px;
  height: 24px;
  display: grid;
}

.msp-viewer-sequence-inspector-dot {
  background: var(--msp-color-text-subtle);
  border-radius: 50%;
  flex: none;
  width: 8px;
  height: 8px;
}

.msp-viewer-sequence-inspector-polymer {
  background: var(--msp-accent-blue);
}

.msp-viewer-sequence-inspector-ligand {
  background: var(--msp-accent-orange);
}

.msp-viewer-sequence-inspector-branched {
  background: var(--msp-accent-green);
}

.msp-viewer-sequence-inspector-text {
  min-width: 0;
  padding: var(--msp-space-1) 0;
  text-align: left;
  flex: 1;
  display: grid;
}

.msp-viewer-sequence-inspector-name, .msp-viewer-sequence-inspector-meta {
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-sequence-inspector-name {
  font-weight: 500;
}

.msp-viewer-sequence-inspector-meta {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-sequence-inspector-chevron {
  color: var(--msp-color-text-subtle);
  flex: none;
  transition: transform .12s ease-out;
  display: grid;
}

.msp-viewer-sequence-inspector-chevronOpen {
  transform: rotate(90deg);
}

.msp-viewer-sequence-inspector-chips {
  gap: var(--msp-space-1);
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-2) var(--msp-space-6);
  flex-wrap: wrap;
  display: flex;
}

.msp-viewer-sequence-inspector-chip {
  min-width: 28px;
  height: 24px;
  padding: 0 var(--msp-space-2);
  border: 1px solid var(--msp-color-border);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-pressed="true"] {
    border-color: var(--msp-color-accent);
    color: var(--msp-color-text);
  }
}

.msp-viewer-sequence-inspector-more {
  align-items: baseline;
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) var(--msp-space-2);
  margin: 0 calc(-1 * var(--msp-space-2));
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-accent);
  font-size: var(--msp-font-size-sm);
  text-align: left;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-sequence-inspector-moreCount {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}
}
`);
