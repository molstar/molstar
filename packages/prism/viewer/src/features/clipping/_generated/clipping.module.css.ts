/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from clipping.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'caption': 'msp-viewer-clipping-caption',
  'cutIcon': 'msp-viewer-clipping-cutIcon',
  'cutRow': 'msp-viewer-clipping-cutRow',
  'cutText': 'msp-viewer-clipping-cutText',
  'fade': 'msp-viewer-clipping-fade',
  'fold': 'msp-viewer-clipping-fold',
  'foldChevron': 'msp-viewer-clipping-foldChevron',
  'foldTitle': 'msp-viewer-clipping-foldTitle',
} as const;

defineStyles('viewer/features/clipping/clipping.module.css', `@layer msp.components {
.msp-viewer-clipping-cutRow {
  align-items: center;
  gap: var(--msp-space-2);
  min-width: 0;
  display: flex;
}

.msp-viewer-clipping-cutIcon {
  margin-right: var(--msp-space-1);
  color: var(--msp-color-text-muted);
  flex: none;
  display: flex;
}

.msp-viewer-clipping-cutText {
  flex: 1;
  min-width: 0;
  display: grid;
}

.msp-viewer-clipping-fade {
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-clipping-caption {
  font-size: var(--msp-font-size-xs);
  color: var(--msp-color-text-muted);
  margin: 0;
}

.msp-viewer-clipping-fold {
  align-items: center;
  gap: var(--msp-space-2);
  width: 100%;
  min-width: 0;
  color: var(--msp-color-text-muted);
  font: inherit;
  font-size: var(--msp-font-size-xs);
  text-align: start;
  cursor: pointer;
  background: none;
  border: none;
  padding: 2px 0;
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-viewer-clipping-foldTitle {
  flex: 1;
  min-width: 0;
}

.msp-viewer-clipping-foldChevron {
  transition: transform var(--msp-motion-fast, .12s) ease;
  flex: none;

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}
}
`);
