/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from measurements-add.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'action': 'msp-viewer-measurements-add-action',
  'add': 'msp-viewer-measurements-add-add',
  'added': 'msp-viewer-measurements-add-added',
  'addedAction': 'msp-viewer-measurements-add-addedAction',
  'addedIcon': 'msp-viewer-measurements-add-addedIcon',
  'addedItem': 'msp-viewer-measurements-add-addedItem',
  'addedRow': 'msp-viewer-measurements-add-addedRow',
  'addedValue': 'msp-viewer-measurements-add-addedValue',
  'caption': 'msp-viewer-measurements-add-caption',
  'chevron': 'msp-viewer-measurements-add-chevron',
  'grip': 'msp-viewer-measurements-add-grip',
  'hint': 'msp-viewer-measurements-add-hint',
  'index': 'msp-viewer-measurements-add-index',
  'kind': 'msp-viewer-measurements-add-kind',
  'kindHint': 'msp-viewer-measurements-add-kindHint',
  'kindTitle': 'msp-viewer-measurements-add-kindTitle',
  'kinds': 'msp-viewer-measurements-add-kinds',
  'more': 'msp-viewer-measurements-add-more',
  'optionsToggle': 'msp-viewer-measurements-add-optionsToggle',
  'point': 'msp-viewer-measurements-add-point',
  'pointActions': 'msp-viewer-measurements-add-pointActions',
  'pointLabel': 'msp-viewer-measurements-add-pointLabel',
  'points': 'msp-viewer-measurements-add-points',
  'slot': 'msp-viewer-measurements-add-slot',
  'style': 'msp-viewer-measurements-add-style',
} as const;

defineStyles('viewer/views/structure/measurements-add.module.css', `@layer msp.components {
.msp-viewer-measurements-add-add {
  font-size: var(--msp-font-size-sm);
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-measurements-add-caption {
  margin: var(--msp-space-1) 0 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-weight: 500;
}

.msp-viewer-measurements-add-points {
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-measurements-add-hint, .msp-viewer-measurements-add-more {
  padding: var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  line-height: 1.45;
}

.msp-viewer-measurements-add-more {
  padding: var(--msp-space-1) var(--msp-space-2);
}

.msp-viewer-measurements-add-point {
  align-items: center;
  gap: var(--msp-space-2);
  height: 26px;
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-measurements-add-slot {
  flex: none;
  place-items: center start;
  min-width: 1.5ch;
  display: grid;

  & > * {
    grid-area: 1 / 1;
  }
}

.msp-viewer-measurements-add-grip {
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

  .msp-viewer-measurements-add-point:hover &, [data-dragging] & {
    opacity: 1;
  }
}

.msp-viewer-measurements-add-point:hover .msp-viewer-measurements-add-index, .msp-viewer-measurements-add-point[data-dragging] .msp-viewer-measurements-add-index {
  visibility: hidden;
}

.msp-viewer-measurements-add-point[data-dragging] {
  opacity: .5;
}

.msp-viewer-measurements-add-point[data-drop="before"] {
  box-shadow: 0 -2px 0 var(--msp-color-text-muted);
}

.msp-viewer-measurements-add-point[data-drop="after"] {
  box-shadow: 0 2px 0 var(--msp-color-text-muted);
}

.msp-viewer-measurements-add-index {
  min-width: 1.5ch;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  flex: none;
}

.msp-viewer-measurements-add-pointLabel {
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

.msp-viewer-measurements-add-pointActions {
  flex: none;
  gap: 1px;
  display: flex;
}

.msp-viewer-measurements-add-action {
  border-radius: var(--msp-radius-sm);
  width: 20px;
  height: 20px;
  color: var(--msp-color-text-muted);
  opacity: 0;
  place-items: center;
  display: grid;

  .msp-viewer-measurements-add-point:hover &, &:focus-visible {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-viewer-measurements-add-kinds {
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-measurements-add-kind {
  align-items: center;
  gap: var(--msp-space-3);
  height: 28px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  text-align: left;
  display: flex;

  &:hover:not(:disabled) {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &:disabled {
    color: var(--msp-color-text-subtle);
    cursor: default;

    & .msp-viewer-measurements-add-kindTitle {
      color: var(--msp-color-text-subtle);
    }
  }
}

.msp-viewer-measurements-add-kindTitle {
  color: var(--msp-color-text);
  flex: none;
  font-weight: 500;
}

.msp-viewer-measurements-add-kindHint {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-align: right;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-measurements-add-optionsToggle {
  align-items: center;
  gap: var(--msp-space-1);
  height: 28px;
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  font-weight: 600;
  display: flex;
}

.msp-viewer-measurements-add-chevron {
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}

.msp-viewer-measurements-add-added {
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-measurements-add-addedItem {
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-measurements-add-addedRow {
  align-items: center;
  gap: var(--msp-space-2);
  height: 26px;
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[data-hidden] .msp-viewer-measurements-add-pointLabel, &[data-hidden] .msp-viewer-measurements-add-addedIcon {
    color: var(--msp-color-text-subtle);
  }
}

.msp-viewer-measurements-add-addedIcon {
  color: var(--msp-color-text-muted);
  flex: none;
}

.msp-viewer-measurements-add-addedValue {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  font-variant-numeric: tabular-nums;
  flex: none;
}

.msp-viewer-measurements-add-addedAction {
  border-radius: var(--msp-radius-sm);
  width: 20px;
  height: 20px;
  color: var(--msp-color-text-muted);
  opacity: 0;
  place-items: center;
  display: grid;

  .msp-viewer-measurements-add-addedRow:hover &, [data-hidden] > .msp-viewer-measurements-add-pointActions > &:first-child, &[aria-expanded="true"], &:focus-visible {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-viewer-measurements-add-style {
  margin: var(--msp-space-1) 0 var(--msp-space-1) calc(var(--msp-space-2) + var(--msp-icon-md) + var(--msp-space-2));
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-1) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  box-shadow: inset 0 0 0 1px var(--msp-glass-edge);
}
}
`);
