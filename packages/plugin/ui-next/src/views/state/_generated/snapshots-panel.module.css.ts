/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from snapshots-panel.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'action': 'msp-state-snapshots-panel-action',
  'actionLabel': 'msp-state-snapshots-panel-actionLabel',
  'actions': 'msp-state-snapshots-panel-actions',
  'block': 'msp-state-snapshots-panel-block',
  'blockCaption': 'msp-state-snapshots-panel-blockCaption',
  'choice': 'msp-state-snapshots-panel-choice',
  'choiceHint': 'msp-state-snapshots-panel-choiceHint',
  'choiceTitle': 'msp-state-snapshots-panel-choiceTitle',
  'description': 'msp-state-snapshots-panel-description',
  'descriptionField': 'msp-state-snapshots-panel-descriptionField',
  'empty': 'msp-state-snapshots-panel-empty',
  'field': 'msp-state-snapshots-panel-field',
  'fields': 'msp-state-snapshots-panel-fields',
  'grip': 'msp-state-snapshots-panel-grip',
  'index': 'msp-state-snapshots-panel-index',
  'inspection': 'msp-state-snapshots-panel-inspection',
  'inspector': 'msp-state-snapshots-panel-inspector',
  'inspectorEmpty': 'msp-state-snapshots-panel-inspectorEmpty',
  'inspectorHeader': 'msp-state-snapshots-panel-inspectorHeader',
  'inspectorKind': 'msp-state-snapshots-panel-inspectorKind',
  'inspectorTitle': 'msp-state-snapshots-panel-inspectorTitle',
  'item': 'msp-state-snapshots-panel-item',
  'key': 'msp-state-snapshots-panel-key',
  'label': 'msp-state-snapshots-panel-label',
  'list': 'msp-state-snapshots-panel-list',
  'menu': 'msp-state-snapshots-panel-menu',
  'menuItem': 'msp-state-snapshots-panel-menuItem',
  'name': 'msp-state-snapshots-panel-name',
  'optionsButton': 'msp-state-snapshots-panel-optionsButton',
  'panel': 'msp-state-snapshots-panel-panel',
  'row': 'msp-state-snapshots-panel-row',
  'rowActions': 'msp-state-snapshots-panel-rowActions',
  'slot': 'msp-state-snapshots-panel-slot',
  'textarea': 'msp-state-snapshots-panel-textarea',
  'thumb': 'msp-state-snapshots-panel-thumb',
} as const;

defineStyles('views/state/snapshots-panel.module.css', `@layer msp.components {
.msp-state-snapshots-panel-panel {
  gap: var(--msp-space-2);
  height: 100%;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  display: flex;
  container-type: inline-size;
}

.msp-state-snapshots-panel-actions {
  gap: var(--msp-space-1);
  flex: none;
  min-width: 0;
  display: flex;
}

@container (width < 280px) {
  .msp-state-snapshots-panel-actions > button:has(.msp-state-snapshots-panel-actionLabel) {
    width: 24px;
    padding: 0;
  }

  .msp-state-snapshots-panel-actionLabel {
    clip-path: inset(50%);
    white-space: nowrap;
    width: 1px;
    height: 1px;
    position: absolute;
    overflow: hidden;
  }
}

.msp-state-snapshots-panel-optionsButton {
  margin-left: auto;
  margin-right: var(--msp-space-2);
}

.msp-state-snapshots-panel-block {
  max-height: 50%;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 50%);
  flex: none;
  display: grid;
  overflow: auto;
}

.msp-state-snapshots-panel-blockCaption {
  padding: var(--msp-space-1) 0 var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-state-snapshots-panel-choice {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  height: 30px;
  margin: 0 calc(-1 * var(--msp-space-2));
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  text-align: left;
  display: flex;

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
    outline: none;
  }
}

.msp-state-snapshots-panel-choiceTitle {
  color: var(--msp-color-text);
  flex: none;
  font-weight: 600;
}

.msp-state-snapshots-panel-choiceHint {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  overflow: hidden;
}

.msp-state-snapshots-panel-list {
  flex: 1 1 0;
  grid-template-columns: minmax(0, 1fr);
  align-content: start;
  min-height: 80px;
  display: grid;
  overflow: auto;
}

.msp-state-snapshots-panel-empty {
  padding: var(--msp-space-3) var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  line-height: 1.5;

  & strong {
    color: var(--msp-color-text-muted);
    font-weight: 600;
  }
}

.msp-state-snapshots-panel-item {
  min-width: 0;
  display: grid;
}

.msp-state-snapshots-panel-row {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  height: 28px;
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  cursor: default;
  user-select: none;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: -2px;
  }

  &[aria-selected="true"] {
    background: var(--msp-color-control);

    & .msp-state-snapshots-panel-name {
      color: var(--msp-color-accent);
    }
  }
}

.msp-state-snapshots-panel-slot {
  flex: none;
  place-items: center end;
  min-width: 2ch;
  display: grid;

  & > * {
    grid-area: 1 / 1;
  }
}

.msp-state-snapshots-panel-grip {
  width: 14px;
  height: 22px;
  color: var(--msp-color-text-subtle);
  cursor: grab;
  touch-action: none;
  opacity: 0;
  place-items: center;
  display: grid;

  &:hover {
    color: var(--msp-color-text);
  }

  .msp-state-snapshots-panel-row:hover &, [data-dragging] & {
    opacity: 1;
  }
}

.msp-state-snapshots-panel-row:hover .msp-state-snapshots-panel-index, [data-dragging] .msp-state-snapshots-panel-index {
  visibility: hidden;
}

.msp-state-snapshots-panel-item {
  &[data-dragging] {
    opacity: .5;
  }

  &[data-drop="before"] {
    box-shadow: 0 -2px 0 var(--msp-color-text-muted);
  }

  &[data-drop="after"] {
    box-shadow: 0 2px 0 var(--msp-color-text-muted);
  }
}

.msp-state-snapshots-panel-index {
  min-width: 2ch;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  text-align: right;
  flex: none;
}

.msp-state-snapshots-panel-thumb {
  border-radius: var(--msp-radius-sm);
  width: 22px;
  height: 16px;
  box-shadow: inset 0 0 0 1px var(--msp-glass-edge);
  flex: none;
}

.msp-state-snapshots-panel-name {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 0 auto;
  min-width: 3ch;
  max-width: 60%;
  font-weight: 500;
  overflow: hidden;
}

.msp-state-snapshots-panel-key {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
}

.msp-state-snapshots-panel-description {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  overflow: hidden;
}

.msp-state-snapshots-panel-rowActions {
  flex: none;
  gap: 1px;
  display: flex;
}

.msp-state-snapshots-panel-action {
  border-radius: var(--msp-radius-sm);
  width: 22px;
  height: 22px;
  color: var(--msp-color-text-muted);
  opacity: 0;
  place-items: center;
  display: grid;

  .msp-state-snapshots-panel-row:hover &, .msp-state-snapshots-panel-row:focus-within &, &[aria-expanded="true"] {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-state-snapshots-panel-menu {
  margin: 2px var(--msp-space-1) var(--msp-space-1) calc(var(--msp-space-2) + 2ch + var(--msp-space-3));
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  display: grid;
}

.msp-state-snapshots-panel-menuItem {
  align-items: center;
  gap: var(--msp-space-3);
  height: 26px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-muted);
  text-align: left;
  grid-template-columns: 14px 1fr;
  display: grid;

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[data-danger] {
    margin-top: var(--msp-space-1);
    border-radius: 0 0 var(--msp-radius-sm) var(--msp-radius-sm);
    box-shadow: 0 -1px 0 var(--msp-glass-edge);
    color: var(--msp-color-danger);

    &:hover, &:focus-visible {
      background: color-mix(in oklab, var(--msp-color-danger), transparent 88%);
      color: var(--msp-color-danger);
    }
  }
}

.msp-state-snapshots-panel-inspector {
  height: clamp(200px, 45%, 440px);
  margin: 0 calc(-1 * var(--msp-space-2));
  border-top: 1px solid var(--msp-glass-edge);
  flex: none;
  overflow: auto;
}

.msp-state-snapshots-panel-inspectorEmpty {
  height: 100%;
  padding: var(--msp-space-4);
  color: var(--msp-color-text-subtle);
  text-align: center;
  place-items: center;
  margin: 0;
  display: grid;
}

.msp-state-snapshots-panel-inspection {
  align-content: start;
  gap: var(--msp-space-3);
  padding: var(--msp-space-3) var(--msp-space-4);
  display: grid;
}

.msp-state-snapshots-panel-inspectorHeader {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 26px;
  display: flex;
}

.msp-state-snapshots-panel-inspectorTitle {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  min-width: 0;
  font-weight: 600;
  overflow: hidden;
}

.msp-state-snapshots-panel-inspectorKind {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-state-snapshots-panel-fields {
  gap: 2px;
  display: grid;
}

.msp-state-snapshots-panel-field {
  grid-template-columns: var(--msp-label-width) minmax(0, 1fr);
  align-items: center;
  gap: var(--msp-space-3);
  min-height: 24px;
  display: grid;
}

.msp-state-snapshots-panel-label {
  color: var(--msp-color-text-muted);
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
}

.msp-state-snapshots-panel-descriptionField {
  gap: var(--msp-space-1);
  padding-top: var(--msp-space-1);
  display: grid;
}

.msp-state-snapshots-panel-textarea {
  min-height: 72px;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font: inherit;
  resize: vertical;
  border: 0;
  line-height: 1.45;

  &::placeholder {
    color: var(--msp-color-text-subtle);
  }

  &:hover {
    background: var(--msp-color-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 1px;
  }
}
}
`);
