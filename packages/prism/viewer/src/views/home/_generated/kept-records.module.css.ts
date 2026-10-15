/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from kept-records.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'asking': 'msp-viewer-kept-records-asking',
  'caption': 'msp-viewer-kept-records-caption',
  'forget': 'msp-viewer-kept-records-forget',
  'group': 'msp-viewer-kept-records-group',
  'listContent': 'msp-viewer-kept-records-listContent',
  'open': 'msp-viewer-kept-records-open',
  'record': 'msp-viewer-kept-records-record',
  'title': 'msp-viewer-kept-records-title',
  'when': 'msp-viewer-kept-records-when',
} as const;

defineStyles('viewer/views/home/kept-records.module.css', `@layer msp.components {
.msp-viewer-kept-records-group {
  gap: var(--msp-space-1);
  display: grid;
}

.msp-viewer-kept-records-caption {
  padding: 0 var(--msp-space-3) var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  font-weight: 500;
}

.msp-viewer-kept-records-listContent {
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-kept-records-record {
  border-radius: var(--msp-radius-md);
  grid-template-columns: minmax(0, 1fr) auto;
  align-items: center;
  display: grid;

  &:hover, &:focus-within {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-kept-records-open {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  height: 26px;
  padding: 0 var(--msp-space-2) 0 var(--msp-inset);
  color: var(--msp-color-text);
  text-align: left;
  grid-template-columns: minmax(0, 1fr) auto;
  display: grid;
}

.msp-viewer-kept-records-title {
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-kept-records-when {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}

.msp-viewer-kept-records-forget {
  visibility: hidden;
  margin-right: var(--msp-space-1);

  .msp-viewer-kept-records-record:is(:hover, :focus-within) & {
    visibility: visible;
  }
}

.msp-viewer-kept-records-asking {
  visibility: visible;
}
}
`);
