/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from contents-summary.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'branched': 'msp-viewer-contents-summary-branched',
  'chevron': 'msp-viewer-contents-summary-chevron',
  'contents': 'msp-viewer-contents-summary-contents',
  'count': 'msp-viewer-contents-summary-count',
  'counts': 'msp-viewer-contents-summary-counts',
  'dot': 'msp-viewer-contents-summary-dot',
  'entities': 'msp-viewer-contents-summary-entities',
  'entity': 'msp-viewer-contents-summary-entity',
  'label': 'msp-viewer-contents-summary-label',
  'ligand': 'msp-viewer-contents-summary-ligand',
  'list': 'msp-viewer-contents-summary-list',
  'meta': 'msp-viewer-contents-summary-meta',
  'more': 'msp-viewer-contents-summary-more',
  'moreCount': 'msp-viewer-contents-summary-moreCount',
  'name': 'msp-viewer-contents-summary-name',
  'other': 'msp-viewer-contents-summary-other',
  'polymer': 'msp-viewer-contents-summary-polymer',
  'text': 'msp-viewer-contents-summary-text',
  'value': 'msp-viewer-contents-summary-value',
  'water': 'msp-viewer-contents-summary-water',
} as const;

defineStyles('viewer/views/structure/contents-summary.module.css', `@layer msp.components {
.msp-viewer-contents-summary-contents {
  gap: var(--msp-space-3);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-contents-summary-counts {
  gap: var(--msp-space-1);
  grid-template-columns: repeat(auto-fit, minmax(64px, 1fr));
  margin: 0;
  display: grid;
}

.msp-viewer-contents-summary-count {
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  flex-direction: column-reverse;
  gap: 2px;
  display: flex;
}

.msp-viewer-contents-summary-value {
  font-size: var(--msp-font-size-lg);
  font-variant-numeric: tabular-nums;
  margin: 0;
  font-weight: 600;
  line-height: 1.1;
}

.msp-viewer-contents-summary-label {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-contents-summary-entities {
  margin: 0 calc(-1 * var(--msp-space-2));
  display: grid;
}

.msp-viewer-contents-summary-list {
  display: grid;
}

.msp-viewer-contents-summary-entity {
  align-items: center;
  gap: var(--msp-space-3);
  min-height: 40px;
  padding: var(--msp-space-1) var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  text-align: left;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-contents-summary-more {
  align-items: baseline;
  gap: var(--msp-space-2);
  padding: var(--msp-space-1) var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-accent);
  font-size: var(--msp-font-size-sm);
  text-align: left;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-contents-summary-moreCount {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-contents-summary-dot {
  background: var(--msp-color-text-subtle);
  border-radius: 50%;
  flex: none;
  width: 10px;
  height: 10px;
}

.msp-viewer-contents-summary-polymer {
  background: var(--msp-accent-blue);
}

.msp-viewer-contents-summary-ligand {
  background: var(--msp-accent-orange);
}

.msp-viewer-contents-summary-branched {
  background: var(--msp-accent-green);
}

.msp-viewer-contents-summary-water {
  background: var(--msp-accent-cyan);
}

.msp-viewer-contents-summary-other {
  background: var(--msp-accent-gray);
}

.msp-viewer-contents-summary-text {
  flex: 1;
  min-width: 0;
  display: grid;
}

.msp-viewer-contents-summary-name, .msp-viewer-contents-summary-meta {
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-contents-summary-name {
  font-weight: 500;
}

.msp-viewer-contents-summary-meta {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-contents-summary-chevron {
  color: var(--msp-color-text-subtle);
  flex: none;
}
}
`);
