/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sequence-chains.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'box': 'msp-viewer-sequence-chains-box',
  'branched': 'msp-viewer-sequence-chains-branched',
  'chains': 'msp-viewer-sequence-chains-chains',
  'chevron': 'msp-viewer-sequence-chains-chevron',
  'chevronOpen': 'msp-viewer-sequence-chains-chevronOpen',
  'chip': 'msp-viewer-sequence-chains-chip',
  'chips': 'msp-viewer-sequence-chains-chips',
  'dot': 'msp-viewer-sequence-chains-dot',
  'entity': 'msp-viewer-sequence-chains-entity',
  'ligand': 'msp-viewer-sequence-chains-ligand',
  'list': 'msp-viewer-sequence-chains-list',
  'meta': 'msp-viewer-sequence-chains-meta',
  'more': 'msp-viewer-sequence-chains-more',
  'moreCount': 'msp-viewer-sequence-chains-moreCount',
  'name': 'msp-viewer-sequence-chains-name',
  'polymer': 'msp-viewer-sequence-chains-polymer',
  'row': 'msp-viewer-sequence-chains-row',
  'text': 'msp-viewer-sequence-chains-text',
} as const;

defineStyles('viewer/views/structure/sequence-chains.module.css', `@layer msp.components {
.msp-viewer-sequence-chains-chains {
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-sequence-chains-list {
  margin: 0 calc(-1 * var(--msp-space-2));
  display: grid;
}

.msp-viewer-sequence-chains-entity {
  min-width: 0;
  display: grid;
}

.msp-viewer-sequence-chains-row {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 36px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-sequence-chains-box {
  border: 1px solid var(--msp-color-border-strong);
  border-radius: var(--msp-radius-sm);
  width: 16px;
  height: 16px;
  color: var(--msp-color-on-accent);
  flex: none;
  place-items: center;
  display: grid;

  &[aria-checked="true"], &[aria-checked="mixed"] {
    border-color: var(--msp-color-accent);
    background: var(--msp-color-accent);
  }
}

.msp-viewer-sequence-chains-dot {
  background: var(--msp-color-text-subtle);
  border-radius: 50%;
  flex: none;
  width: 8px;
  height: 8px;
}

.msp-viewer-sequence-chains-polymer {
  background: var(--msp-accent-blue);
}

.msp-viewer-sequence-chains-ligand {
  background: var(--msp-accent-orange);
}

.msp-viewer-sequence-chains-branched {
  background: var(--msp-accent-green);
}

.msp-viewer-sequence-chains-text {
  min-width: 0;
  padding: var(--msp-space-1) 0;
  text-align: left;
  flex: 1;
  display: grid;
}

.msp-viewer-sequence-chains-name, .msp-viewer-sequence-chains-meta {
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-sequence-chains-name {
  font-weight: 500;
}

.msp-viewer-sequence-chains-meta {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-sequence-chains-chevron {
  color: var(--msp-color-text-subtle);
  flex: none;
  transition: transform .12s ease-out;
  display: grid;
}

.msp-viewer-sequence-chains-chevronOpen {
  transform: rotate(90deg);
}

.msp-viewer-sequence-chains-chips {
  gap: var(--msp-space-1);
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-2) var(--msp-space-6);
  flex-wrap: wrap;
  display: flex;
}

.msp-viewer-sequence-chains-chip {
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

.msp-viewer-sequence-chains-more {
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

.msp-viewer-sequence-chains-moreCount {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}
}
`);
