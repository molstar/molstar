/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from shell.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'busy': 'msui-shell-busy',
  'busyCard': 'msui-shell-busyCard',
  'busyRow': 'msui-shell-busyRow',
  'card': 'msui-shell-card',
  'cardSlot': 'msui-shell-cardSlot',
  'compact': 'msui-shell-compact',
  'dock': 'msui-shell-dock',
  'drop': 'msui-shell-drop',
  'grip': 'msui-shell-grip',
  'group': 'msui-shell-group',
  'main': 'msui-shell-main',
  'mainContent': 'msui-shell-mainContent',
  'modal': 'msui-shell-modal',
  'node': 'msui-shell-node',
  'page': 'msui-shell-page',
  'progress': 'msui-shell-progress',
  'rail': 'msui-shell-rail',
  'railButton': 'msui-shell-railButton',
  'railCaption': 'msui-shell-railCaption',
  'railFoot': 'msui-shell-railFoot',
  'railSection': 'msui-shell-railSection',
  'railToggle': 'msui-shell-railToggle',
  'scrim': 'msui-shell-scrim',
  'sheet': 'msui-shell-sheet',
  'sheetBar': 'msui-shell-sheetBar',
  'sheetButton': 'msui-shell-sheetButton',
  'sheetContent': 'msui-shell-sheetContent',
  'sheetTabs': 'msui-shell-sheetTabs',
  'shell': 'msui-shell-shell',
  'split': 'msui-shell-split',
  'splitter': 'msui-shell-splitter',
  'status': 'msui-shell-status',
  'tab': 'msui-shell-tab',
  'tabs': 'msui-shell-tabs',
  'window': 'msui-shell-window',
  'windowTabs': 'msui-shell-windowTabs',
  'windows': 'msui-shell-windows',
} as const;

defineStyles('shell/shell.module.css', `@layer msui.components {
.msui-shell-shell {
  --msui-shell-rail: calc(10 * var(--msui-unit));
  --msui-shell-button: calc(8 * var(--msui-unit));
  --msui-shell-caption: calc(3 * var(--msui-unit));
  --msui-shell-head: calc(var(--msui-surface-pad) + var(--msui-shell-button));
  --msui-shell-splitter: calc(2 * var(--msui-unit));
  --msui-shell-line: color-mix(in oklab, var(--msui-line) 60%, transparent);
  background: var(--msui-bg);
  width: 100%;
  min-width: 0;
  height: 100%;
  min-height: 0;
  color: var(--msui-text);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
  display: flex;
  position: relative;
  overflow: hidden;
  container: shell / size;
}

@container shell (height <= 619px) {
  .msui-shell-rail, .msui-shell-dock {
    --msui-shell-button: calc(7 * var(--msui-unit));
  }
}

@container shell (height <= 519px) {
  .msui-shell-rail, .msui-shell-dock {
    --msui-shell-button: calc(6 * var(--msui-unit));
  }
}

.msui-shell-dock {
  --msui-shell-head: calc(var(--msui-surface-pad) + var(--msui-shell-button));
  flex: auto;
  min-width: 0;
  min-height: 0;
  display: flex;
}

.msui-shell-split {
  flex: auto;
  min-width: 0;
  min-height: 0;
  display: flex;
  position: relative;

  &[data-dir="row"] {
    flex-direction: row;
  }

  &[data-dir="column"] {
    flex-direction: column;
  }
}

.msui-shell-node {
  min-width: 0;
  min-height: 0;
  display: flex;
  position: relative;

  &[data-flex] {
    flex: 1 1 0;
  }

  &:not([data-flex]) {
    flex: 0 0 var(--msui-shell-size);
  }

  &[data-shown="false"] {
    display: none;
  }

  &[data-folded] {
    flex: none;

    & [data-msui-view] > :not(header) {
      display: none;
    }
  }
}

.msui-shell-main {
  flex: auto;
  min-width: 0;
  min-height: 0;
  position: relative;
  overflow: hidden;
  container: shell-main / size;
}

.msui-shell-mainContent {
  display: flex;
  position: absolute;
  inset: 0;
}

.msui-shell-status {
  bottom: var(--msui-space-2);
  left: var(--msui-space-2);
  max-width: calc(100% - 2 * var(--msui-space-2));
  pointer-events: none;
  position: absolute;

  & > * {
    pointer-events: auto;
  }
}

.msui-shell-railToggle {
  top: var(--msui-surface-pad);
  right: var(--msui-surface-pad);
  z-index: 1;
  position: absolute;
}

.msui-shell-group {
  background: var(--msui-surface);
  min-width: 0;
  min-height: 0;
  box-shadow: inset 1px 0 0 var(--msui-shell-line);
  flex-direction: column;
  flex: auto;
  display: flex;

  [data-dir="column"] > .msui-shell-node + .msui-shell-node > & {
    box-shadow: inset 1px 0 0 var(--msui-shell-line),
      inset 0 1px 0 var(--msui-shell-line);
  }
}

.msui-shell-splitter {
  z-index: 2;
  touch-action: none;
  background: none;
  position: absolute;

  &[data-dir="row"] {
    width: var(--msui-shell-splitter);
    cursor: col-resize;
    top: 0;
    bottom: 0;
  }

  &[data-dir="column"] {
    height: var(--msui-shell-splitter);
    cursor: row-resize;
    left: 0;
    right: 0;
  }

  &[data-at="start"][data-dir="row"] {
    left: calc(-.5 * var(--msui-shell-splitter));
  }

  &[data-at="end"][data-dir="row"] {
    right: calc(-.5 * var(--msui-shell-splitter));
  }

  &[data-at="start"][data-dir="column"] {
    top: calc(-.5 * var(--msui-shell-splitter));
  }

  &[data-at="end"][data-dir="column"] {
    bottom: calc(-.5 * var(--msui-shell-splitter));
  }

  &:hover, &:focus-visible, &[data-dragging] {
    background: color-mix(in oklab, var(--msui-accent) 40%, transparent);
    outline: none;
  }
}

.msui-shell-rail {
  width: var(--msui-shell-rail);
  min-height: 0;
  padding: var(--msui-surface-pad);
  background: var(--msui-surface);
  box-shadow: inset 1px 0 0 var(--msui-shell-line);
  flex-direction: column;
  flex: none;
  align-items: stretch;
  display: flex;
}

.msui-shell-railSection {
  flex-direction: column;
  gap: 0;
  display: flex;
  position: relative;

  & + & {
    margin-top: var(--msui-space-1);
    padding-top: var(--msui-space-1);

    &:before {
      content: "";
      top: 0;
      left: var(--msui-space-1);
      right: var(--msui-space-1);
      background: var(--msui-shell-line);
      height: 1px;
      position: absolute;
    }
  }
}

.msui-shell-railFoot > .msui-shell-railSection + .msui-shell-railSection {
  margin-top: 0;
  padding-top: 0;

  &:before {
    content: none;
  }
}

.msui-shell-railFoot {
  padding-top: var(--msui-space-1);
  flex-direction: column;
  margin-top: auto;
  display: flex;
}

.msui-shell-railButton {
  --msui-icon-size: var(--msui-icon);
  width: 100%;
  height: var(--msui-shell-button);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  transition: background-color var(--msui-motion-fast) var(--msui-ease),
    color var(--msui-motion-fast) var(--msui-ease);
  flex-direction: column;
  justify-content: center;
  align-items: center;
  display: flex;

  &:hover {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &[aria-selected="true"], &[aria-pressed="true"] {
    background: var(--msui-control-active);
    color: var(--msui-text);
  }

  &[data-captioned] {
    height: calc(var(--msui-shell-button) + var(--msui-shell-caption));
  }

  &[hidden] {
    display: none;
  }
}

.msui-shell-railCaption {
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-shell-caption);
}

.msui-shell-windows {
  z-index: 2;
  pointer-events: none;
  position: absolute;
  inset: 0;
}

.msui-shell-window {
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface);
  min-width: 0;
  min-height: 0;
  box-shadow: var(--msui-edge), var(--msui-shadow-float);
  pointer-events: auto;
  flex-direction: column;
  display: flex;
  position: absolute;
  overflow: hidden;

  & header {
    cursor: grab;
    touch-action: none;
  }

  &[data-enter] {
    --msui-rise-y: var(--msui-rise);
    animation: msui-float-in var(--msui-enter);
  }

  &[data-adjusting] {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-shell-windowTabs, .msui-shell-tabs {
  align-items: stretch;
  gap: var(--msui-space-1);
  min-width: 0;
  height: 100%;
  scrollbar-width: none;
  flex: 1;
  margin-inline-start: calc(-1 * var(--msui-inner));
  display: flex;
  overflow-x: auto;
}

.msui-shell-tab {
  align-items: center;
  gap: var(--msui-space-1);
  padding-inline: var(--msui-space-2);
  color: var(--msui-text-muted);
  font-weight: var(--msui-weight-strong);
  white-space: nowrap;
  cursor: pointer;
  --msui-icon-size: var(--msui-icon-sm);
  flex: none;
  display: flex;

  &:hover {
    color: var(--msui-text);
  }

  &[aria-selected="true"] {
    color: var(--msui-text);
    box-shadow: inset 0 -2px 0 var(--msui-accent);
  }
}

.msui-shell-grip {
  z-index: 1;
  touch-action: none;
  position: absolute;

  &[data-edge="right"] {
    top: 0;
    right: 0;
    bottom: var(--msui-space-3);
    width: var(--msui-space-1);
    cursor: ew-resize;
  }

  &[data-edge="bottom"] {
    left: 0;
    right: var(--msui-space-3);
    height: var(--msui-space-1);
    cursor: ns-resize;
    bottom: 0;
  }

  &[data-edge="corner"] {
    width: var(--msui-space-3);
    height: var(--msui-space-3);
    cursor: nwse-resize;
    bottom: 0;
    right: 0;
  }
}

.msui-shell-cardSlot {
  display: contents;
}

.msui-shell-card {
  z-index: 3;
  width: var(--msui-shell-card, 340px);
  background: var(--msui-surface);
  max-width: 100%;
  box-shadow: inset 1px 0 0 var(--msui-shell-line);
  display: flex;
  position: absolute;
  top: 0;
  bottom: 0;
  right: 0;
}

.msui-shell-page {
  z-index: 4;
  background: var(--msui-surface);
  display: flex;
  position: absolute;
  inset: 0;
}

.msui-shell-scrim {
  z-index: var(--msui-z-dialog);
  padding: var(--msui-space-6);
  background: var(--msui-scrim);
  place-items: center;
  display: grid;
  position: absolute;
  inset: 0;
}

.msui-shell-modal {
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface);
  width: min(480px, 100%);
  max-height: 100%;
  box-shadow: var(--msui-edge), var(--msui-shadow-overlay);
  display: flex;
  overflow: hidden;
}

.msui-shell-busy {
  z-index: var(--msui-z-dialog);
  background: var(--msui-scrim);
  place-items: center;
  display: grid;
  position: absolute;
  inset: 0;
}

.msui-shell-busyCard {
  gap: var(--msui-space-2);
  width: min(320px, calc(100% - 2 * var(--msui-space-4)));
  padding: var(--msui-space-3) var(--msui-inset);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: var(--msui-edge), var(--msui-shadow-overlay);
  display: grid;
}

.msui-shell-busyRow {
  align-items: center;
  gap: var(--msui-space-2);
  display: flex;
}

.msui-shell-progress {
  height: var(--msui-space-1);
  border-radius: var(--msui-radius-sm);
  background: var(--msui-control);
  overflow: hidden;

  & > span {
    background: var(--msui-accent);
    height: 100%;
    transition: width var(--msui-motion-fast) var(--msui-ease);
    display: block;
  }

  &[data-indeterminate] > span {
    animation: msui-shell-indeterminate 1.2s var(--msui-ease-in-out) infinite;
    width: 30% !important;
  }
}

@keyframes msui-shell-indeterminate {
  from {
    translate: -100%;
  }

  to {
    translate: 340%;
  }
}

.msui-shell-drop {
  inset: var(--msui-space-2);
  z-index: var(--msui-z-overlay);
  border: 2px dashed var(--msui-accent);
  border-radius: var(--msui-radius-lg);
  background: color-mix(in oklab, var(--msui-bg) 80%, transparent);
  color: var(--msui-text);
  font-size: var(--msui-text-md);
  font-weight: var(--msui-weight-strong);
  pointer-events: none;
  place-items: center;
  display: grid;
  position: absolute;
}

.msui-shell-compact {
  flex-direction: column;
  flex: auto;
  min-width: 0;
  min-height: 0;
  display: flex;
}

.msui-shell-sheet {
  background: var(--msui-surface);
  min-height: 0;
  box-shadow: inset 0 1px 0 var(--msui-shell-line);
  flex-direction: column;
  flex: 0 0 min(45cqh, 420px);
  display: flex;

  &[data-open="false"] {
    flex-basis: auto;
  }
}

.msui-shell-sheetContent {
  flex: auto;
  min-height: 0;
  display: flex;

  &[hidden] {
    display: none;
  }
}

.msui-shell-sheetBar {
  --msui-icon-size: var(--msui-icon);
  height: calc(11 * var(--msui-unit));
  padding-inline: var(--msui-surface-pad);
  box-shadow: inset 0 1px 0 var(--msui-shell-line);
  flex: none;
  align-items: stretch;
  display: flex;
}

.msui-shell-sheetTabs {
  display: contents;
}

.msui-shell-sheetButton {
  min-width: 0;
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  flex-direction: column;
  flex: 1 1 0;
  justify-content: center;
  align-items: center;
  gap: 0;
  display: flex;

  &[aria-selected="true"] {
    color: var(--msui-accent);
  }
}
}
`);
