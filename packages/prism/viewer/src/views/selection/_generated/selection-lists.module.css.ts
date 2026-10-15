/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from selection-lists.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'action': 'msp-viewer-selection-lists-action',
  'actions': 'msp-viewer-selection-lists-actions',
  'block': 'msp-viewer-selection-lists-block',
  'caption': 'msp-viewer-selection-lists-caption',
  'clear': 'msp-viewer-selection-lists-clear',
  'count': 'msp-viewer-selection-lists-count',
  'detail': 'msp-viewer-selection-lists-detail',
  'editor': 'msp-viewer-selection-lists-editor',
  'head': 'msp-viewer-selection-lists-head',
  'item': 'msp-viewer-selection-lists-item',
  'label': 'msp-viewer-selection-lists-label',
  'lead': 'msp-viewer-selection-lists-lead',
  'list': 'msp-viewer-selection-lists-list',
  'listContent': 'msp-viewer-selection-lists-listContent',
  'name': 'msp-viewer-selection-lists-name',
  'part': 'msp-viewer-selection-lists-part',
  'partLine': 'msp-viewer-selection-lists-partLine',
  'partName': 'msp-viewer-selection-lists-partName',
  'row': 'msp-viewer-selection-lists-row',
  'summary': 'msp-viewer-selection-lists-summary',
  'summaryAction': 'msp-viewer-selection-lists-summaryAction',
  'summaryHead': 'msp-viewer-selection-lists-summaryHead',
  'summaryHint': 'msp-viewer-selection-lists-summaryHint',
  'summaryIcon': 'msp-viewer-selection-lists-summaryIcon',
  'summaryLevel': 'msp-viewer-selection-lists-summaryLevel',
  'summaryOffer': 'msp-viewer-selection-lists-summaryOffer',
  'summaryStart': 'msp-viewer-selection-lists-summaryStart',
  'summaryTitle': 'msp-viewer-selection-lists-summaryTitle',
  'swatch': 'msp-viewer-selection-lists-swatch',
  'value': 'msp-viewer-selection-lists-value',
} as const;

defineStyles('viewer/views/selection/selection-lists.module.css', `@layer msp.components {
.msp-viewer-selection-lists-block {
  margin-top: var(--msp-space-3);
  font-size: var(--msp-font-size-sm);
  gap: var(--msp-space-1);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-selection-lists-head {
  align-items: center;
  gap: var(--msp-space-2);
  height: 22px;
  display: flex;
}

.msp-viewer-selection-lists-caption {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: 1;
  font-weight: 500;
}

.msp-viewer-selection-lists-count {
  margin-left: var(--msp-space-1);
  font-variant-numeric: tabular-nums;
}

.msp-viewer-selection-lists-clear {
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
  flex: none;
  font-weight: 500;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-viewer-selection-lists-list {
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  contain: content;
  max-height: 224px;
}

.msp-viewer-selection-lists-listContent {
  padding: var(--msp-space-1);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-selection-lists-item {
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-selection-lists-row {
  align-items: center;
  gap: var(--msp-space-2);
  height: 26px;
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-selection-lists-lead {
  width: 14px;
  color: var(--msp-color-text-muted);
  flex: none;
  place-items: center;
  display: grid;
}

.msp-viewer-selection-lists-swatch {
  width: 12px;
  height: 12px;
  box-shadow: inset 0 0 0 1px var(--msp-glass-edge);
  border-radius: 50%;
}

.msp-viewer-selection-lists-label {
  min-width: 0;
  color: var(--msp-color-text);
  text-align: left;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  font-weight: 500;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-selection-lists-value {
  max-width: 50%;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
  flex: none;
  overflow: hidden;
}

.msp-viewer-selection-lists-name {
  max-width: 45%;
  color: var(--msp-color-text);
  white-space: nowrap;
  flex: none;
  font-weight: 500;
  overflow: hidden;
}

.msp-viewer-selection-lists-detail {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-selection-lists-actions {
  flex: none;
  gap: 1px;
  display: flex;
}

.msp-viewer-selection-lists-action {
  border-radius: var(--msp-radius-sm);
  width: 20px;
  height: 20px;
  color: var(--msp-color-text-muted);
  opacity: 0;
  place-items: center;
  display: grid;

  .msp-viewer-selection-lists-row:hover &, &[aria-expanded="true"], &:focus-visible {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-viewer-selection-lists-editor {
  margin: var(--msp-space-1) 0 var(--msp-space-1) var(--msp-space-5);
  padding: var(--msp-space-1) var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
}

.msp-viewer-selection-lists-summary {
  margin: 0 var(--msp-space-2);
  padding: var(--msp-space-2) var(--msp-space-2) var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  gap: var(--msp-space-1);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-selection-lists-summaryHead {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 24px;
  display: flex;
}

.msp-viewer-selection-lists-summaryIcon {
  color: var(--msp-color-text-muted);
  flex: none;
}

.msp-viewer-selection-lists-summaryTitle {
  min-width: 0;
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-md);
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  font-weight: 600;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-selection-lists-summaryAction {
  border-radius: var(--msp-radius-sm);
  width: 24px;
  height: 24px;
  color: var(--msp-color-text-muted);
  flex: none;
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-viewer-selection-lists-part {
  align-items: baseline;
  gap: var(--msp-space-2);
  min-width: 0;
  padding-left: calc(var(--msp-icon-md) + var(--msp-space-2));
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-viewer-selection-lists-partName {
  max-width: 50%;
  color: var(--msp-color-text-muted);
  white-space: nowrap;
  flex: none;
  font-weight: 500;
  overflow: hidden;
}

.msp-viewer-selection-lists-partLine {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-selection-lists-summaryHint {
  padding-left: calc(var(--msp-icon-md) + var(--msp-space-2));
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  line-height: 1.45;
}

.msp-viewer-selection-lists-summaryLevel {
  color: var(--msp-color-text-muted);
  font-weight: 500;
}

.msp-viewer-selection-lists-summaryOffer {
  justify-items: start;
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-selection-lists-summaryStart {
  align-items: center;
  gap: var(--msp-space-1);
  height: 26px;
  margin-left: calc(var(--msp-icon-md) + var(--msp-space-2));
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-accent);
  font-size: var(--msp-font-size-sm);
  font-weight: 600;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: -2px;
  }
}
}
`);
