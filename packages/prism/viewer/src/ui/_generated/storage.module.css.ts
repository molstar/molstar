/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from storage.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'body': 'msp-viewer-storage-body',
  'caption': 'msp-viewer-storage-caption',
  'chevron': 'msp-viewer-storage-chevron',
  'confirm': 'msp-viewer-storage-confirm',
  'confirmText': 'msp-viewer-storage-confirmText',
  'detail': 'msp-viewer-storage-detail',
  'entry': 'msp-viewer-storage-entry',
  'error': 'msp-viewer-storage-error',
  'footer': 'msp-viewer-storage-footer',
  'header': 'msp-viewer-storage-header',
  'label': 'msp-viewer-storage-label',
  'meta': 'msp-viewer-storage-meta',
  'none': 'msp-viewer-storage-none',
  'section': 'msp-viewer-storage-section',
  'size': 'msp-viewer-storage-size',
  'source': 'msp-viewer-storage-source',
  'sourceHead': 'msp-viewer-storage-sourceHead',
  'sourceTitle': 'msp-viewer-storage-sourceTitle',
  'title': 'msp-viewer-storage-title',
  'toggle': 'msp-viewer-storage-toggle',
  'where': 'msp-viewer-storage-where',
} as const;

defineStyles('viewer/ui/storage.module.css', `@layer msp.components {
.msp-viewer-storage-section {
  border-radius: var(--msp-radius-md);

  &:has(.msp-viewer-storage-toggle[aria-expanded="true"]) {
    background: var(--msp-color-well);
  }
}

.msp-viewer-storage-header {
  border-radius: var(--msp-radius-md);
  align-items: center;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-storage-toggle {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  height: 32px;
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-muted);
  text-align: left;
  flex: 1;
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-viewer-storage-title {
  color: var(--msp-color-text);
  flex: 1;
  font-weight: 600;
}

.msp-viewer-storage-chevron {
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="false"] > & {
    transform: rotate(-90deg);
  }
}

.msp-viewer-storage-body {
  gap: var(--msp-space-3);
  padding: var(--msp-space-1) var(--msp-space-3) var(--msp-space-3);
  display: grid;
}

.msp-viewer-storage-where, .msp-viewer-storage-none, .msp-viewer-storage-error {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  overflow-wrap: anywhere;
  margin: 0;
}

.msp-viewer-storage-error {
  color: var(--msp-color-danger);
}

.msp-viewer-storage-source {
  gap: var(--msp-space-1);
  display: grid;
}

.msp-viewer-storage-sourceHead {
  align-items: baseline;
  gap: var(--msp-space-3);
  display: flex;
}

.msp-viewer-storage-sourceTitle {
  font-weight: 600;
}

.msp-viewer-storage-meta {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin-left: auto;
}

.msp-viewer-storage-entry {
  align-items: center;
  gap: 0 var(--msp-space-2);
  padding: var(--msp-space-1) 0 var(--msp-space-1) var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  grid-template-columns: minmax(0, 1fr) auto auto;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[data-stale] .msp-viewer-storage-detail {
    color: var(--msp-color-warning);
  }
}

.msp-viewer-storage-label {
  overflow-wrap: anywhere;
  min-width: 0;
}

.msp-viewer-storage-size {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
}

.msp-viewer-storage-detail {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  overflow-wrap: anywhere;
  grid-column: 1 / -1;
}

.msp-viewer-storage-footer, .msp-viewer-storage-confirm {
  align-items: center;
  gap: var(--msp-space-2);
  flex-wrap: wrap;
  display: flex;
}

.msp-viewer-storage-confirmText {
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-storage-caption {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin-left: auto;
}
}
`);
