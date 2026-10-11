/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from canvas.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'bottom': 'msp-shell-canvas-bottom',
  'bottomGrid': 'msp-shell-canvas-bottomGrid',
  'card': 'msp-shell-canvas-card',
  'card-in': 'msp-shell-canvas-card-in',
  'cardBody': 'msp-shell-canvas-cardBody',
  'cardHeader': 'msp-shell-canvas-cardHeader',
  'cardTitle': 'msp-shell-canvas-cardTitle',
  'chip': 'msp-shell-canvas-chip',
  'context': 'msp-shell-canvas-context',
  'frame': 'msp-shell-canvas-frame',
  'lead': 'msp-shell-canvas-lead',
  'menuMark': 'msp-shell-canvas-menuMark',
  'message': 'msp-shell-canvas-message',
  'pill': 'msp-shell-canvas-pill',
  'pill-in': 'msp-shell-canvas-pill-in',
  'pillButton': 'msp-shell-canvas-pillButton',
  'progress': 'msp-shell-canvas-progress',
  'sections': 'msp-shell-canvas-sections',
  'spin': 'msp-shell-canvas-spin',
  'stack': 'msp-shell-canvas-stack',
  'status': 'msp-shell-canvas-status',
  'text': 'msp-shell-canvas-text',
  'tool': 'msp-shell-canvas-tool',
  'tool-menu-in': 'msp-shell-canvas-tool-menu-in',
  'toolMenu': 'msp-shell-canvas-toolMenu',
  'toolMenuCard': 'msp-shell-canvas-toolMenuCard',
  'toolMenuCheck': 'msp-shell-canvas-toolMenuCheck',
  'toolMenuHint': 'msp-shell-canvas-toolMenuHint',
  'toolMenuLabel': 'msp-shell-canvas-toolMenuLabel',
  'toolMenuRow': 'msp-shell-canvas-toolMenuRow',
  'toolbar': 'msp-shell-canvas-toolbar',
  'toolbarSlot': 'msp-shell-canvas-toolbarSlot',
} as const;

defineStyles('views/shell/canvas.module.css', `@layer msp.components {
.msp-shell-canvas-frame {
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  width: 44px;
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  flex-direction: column;
  display: flex;
  overflow: hidden;

  & > * + *, & .msp-shell-canvas-sections > * + * {
    border-top: 1px solid var(--msp-glass-edge);
  }

  & .msp-shell-canvas-sections:empty {
    display: none;
  }
}

.msp-shell-canvas-sections, .msp-shell-canvas-toolbar {
  flex-direction: column;
  display: flex;
}

.msp-shell-canvas-tool {
  width: 100%;
  height: 40px;
  color: var(--msp-color-text-muted);
  transition: background-color var(--msp-motion-fast) var(--msp-ease),
    color var(--msp-motion-fast) var(--msp-ease);
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[aria-pressed="true"] {
    background: var(--msp-color-control);
    color: var(--msp-color-accent);
  }
}

.msp-shell-canvas-card {
  top: var(--msp-space-4);
  right: var(--msp-canvas-panel-right, 64px);
  bottom: var(--msp-space-4);
  z-index: var(--msp-z-overlay);
  width: min(var(--msp-right-panel-width, 340px), calc(100% - 96px));
  border-radius: var(--msp-radius-lg);
  animation: msp-shell-canvas-card-in var(--msp-motion-base) var(--msp-ease);
  background: var(--msp-glass-panel-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  flex-direction: column;
  display: flex;
  position: absolute;
  overflow: hidden;

  [data-right-panel="open"] > & {
    animation: none;
  }
}

@keyframes msp-shell-canvas-card-in {
  from {
    opacity: 0;
    transform: translateX(8px);
  }
}

.msp-shell-canvas-cardHeader {
  align-items: center;
  gap: var(--msp-space-3);
  height: 40px;
  padding: 0 var(--msp-space-3) 0 var(--msp-space-4);
  font-weight: 600;
  display: flex;
}

.msp-shell-canvas-cardTitle {
  flex: 1;
}

.msp-shell-canvas-cardBody {
  min-height: 0;
  overflow: auto;
}

.msp-shell-canvas-stack {
  --msp-stack-row: 34px;
  top: var(--msp-canvas-top-inset, var(--msp-space-4));
  right: var(--msp-canvas-right-inset, var(--msp-space-4));
  position: absolute;
}

:where([data-mode="compact"]) .msp-shell-canvas-stack {
  right: var(--msp-canvas-rail-inset, var(--msp-space-4));
}

.msp-shell-canvas-stack .msp-shell-canvas-frame {
  width: 36px;
}

.msp-shell-canvas-stack .msp-shell-canvas-tool {
  height: var(--msp-stack-row);
}

.msp-shell-canvas-stack :is(.msp-shell-canvas-tool, button) > svg:first-child {
  width: 18px;
  height: 18px;
}

.msp-shell-canvas-tool {
  position: relative;
}

.msp-shell-canvas-menuMark {
  opacity: .45;
  border-top: 4px solid #0000;
  border-right: 4px solid;
  position: absolute;
  bottom: 5px;
  right: 5px;
}

.msp-shell-canvas-toolMenu {
  z-index: var(--msp-z-overlay);
  padding-right: var(--msp-space-3);
  animation: msp-shell-canvas-tool-menu-in var(--msp-motion-fast) var(--msp-ease);
  position: absolute;
  right: 100%;
}

@keyframes msp-shell-canvas-tool-menu-in {
  from {
    opacity: 0;
    transform: translateX(4px);
  }
}

.msp-shell-canvas-toolMenuCard {
  min-width: 200px;
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-panel-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  display: grid;
}

.msp-shell-canvas-toolMenuRow {
  align-items: center;
  gap: var(--msp-space-3);
  min-height: 32px;
  padding: var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  text-align: left;
  grid-template-columns: 16px 1fr auto;
  display: grid;

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-shell-canvas-toolMenuLabel {
  color: var(--msp-color-text);
  display: grid;
}

.msp-shell-canvas-toolMenuCheck {
  color: var(--msp-color-accent);
}

.msp-shell-canvas-toolMenuHint {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-shell-canvas-toolbarSlot {
  top: var(--msp-space-4);
  right: var(--msp-space-4);
  position: absolute;
}

.msp-shell-canvas-lead {
  top: var(--msp-canvas-top-inset, var(--msp-space-4));
  left: var(--msp-space-4);
  gap: var(--msp-space-3);
  width: min(320px, 100% - 96px);
  display: grid;
  position: absolute;
}

:where([data-mode="compact"]) .msp-shell-canvas-lead {
  width: min(320px, 100% - 180px);
}

.msp-shell-canvas-bottom {
  left: var(--msp-space-4);
  right: var(--msp-canvas-right-inset, var(--msp-space-4));
  bottom: var(--msp-space-4);
  pointer-events: none;
  position: absolute;
  container: bottom-bar / inline-size;
}

.msp-shell-canvas-bottomGrid {
  align-items: end;
  gap: var(--msp-space-3);
  grid-template-columns: minmax(0, 1fr) auto;
  grid-template-areas: "status context";
  display: grid;

  & > * {
    pointer-events: auto;
  }
}

@container bottom-bar (width < 480px) {
  .msp-shell-canvas-bottomGrid {
    grid-template-columns: minmax(0, 1fr);
    grid-template-areas: "context"
                         "status";
  }
}

.msp-shell-canvas-context {
  justify-content: flex-end;
  justify-self: end;
  gap: var(--msp-space-2);
  flex-wrap: wrap;
  grid-area: context;
  display: flex;
}

.msp-shell-canvas-chip {
  height: 30px;
  padding: 0 var(--msp-space-4);
  border-radius: var(--msp-radius-lg);
  font-size: var(--msp-font-size-md);
}

.msp-shell-canvas-status {
  align-items: center;
  gap: var(--msp-space-2);
  pointer-events: none;
  flex-wrap: wrap;
  grid-area: status;
  display: flex;

  & > * {
    pointer-events: auto;
  }
}

.msp-shell-canvas-pill {
  align-items: center;
  gap: var(--msp-space-2);
  max-width: 100%;
  height: 26px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-lg);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  white-space: nowrap;
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  display: inline-flex;
}

.msp-shell-canvas-message {
  animation: msp-shell-canvas-pill-in var(--msp-motion-base) var(--msp-ease);

  &[data-kind="error"] {
    color: var(--msp-color-danger);
  }

  &[data-kind="warning"] {
    color: var(--msp-color-warning);
  }

  &[data-kind="success"] {
    color: var(--msp-color-success);
  }
}

@keyframes msp-shell-canvas-pill-in {
  from {
    opacity: 0;
    transform: translateY(6px);
  }
}

.msp-shell-canvas-text {
  text-overflow: ellipsis;
  overflow: hidden;
}

.msp-shell-canvas-progress {
  background: var(--msp-color-border);
  border-radius: 2px;
  width: 48px;
  height: 4px;
  overflow: hidden;

  & > div {
    background: var(--msp-color-accent);
    height: 100%;
  }
}

.msp-shell-canvas-spin {
  animation: 1s linear infinite msp-shell-canvas-spin;
}

@keyframes msp-shell-canvas-spin {
  to {
    transform: rotate(360deg);
  }
}

.msp-shell-canvas-pillButton {
  margin-right: calc(-1 * var(--msp-space-2));
  color: var(--msp-color-text-muted);
  border-radius: 50%;
  place-items: center;
  display: grid;

  &:hover {
    color: var(--msp-color-text);
  }
}
}
`);
