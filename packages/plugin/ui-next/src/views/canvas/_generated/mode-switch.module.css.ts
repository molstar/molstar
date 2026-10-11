/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from mode-switch.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'anchor': 'msp-canvas-mode-switch-anchor',
  'chevron': 'msp-canvas-mode-switch-chevron',
  'group': 'msp-canvas-mode-switch-group',
  'label': 'msp-canvas-mode-switch-label',
  'level': 'msp-canvas-mode-switch-level',
  'level-in': 'msp-canvas-mode-switch-level-in',
  'levelButton': 'msp-canvas-mode-switch-levelButton',
  'levelCaption': 'msp-canvas-mode-switch-levelCaption',
  'levelRow': 'msp-canvas-mode-switch-levelRow',
  'levels': 'msp-canvas-mode-switch-levels',
  'segment': 'msp-canvas-mode-switch-segment',
  'switch': 'msp-canvas-mode-switch-switch',
} as const;

defineStyles('views/canvas/mode-switch.module.css', `@layer msp.components {
.msp-canvas-mode-switch-anchor {
  top: var(--msp-canvas-top-inset, var(--msp-space-4));
  right: var(--msp-canvas-right-inset, var(--msp-space-4));
  align-items: flex-end;
  gap: var(--msp-space-3);
  transition: right var(--msp-motion-base) var(--msp-ease);
  flex-direction: column;
  display: flex;
  position: absolute;
}

.msp-canvas-mode-switch-switch, .msp-canvas-mode-switch-levelButton, .msp-canvas-mode-switch-levels {
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
}

.msp-canvas-mode-switch-switch {
  align-items: center;
  gap: 2px;
  height: 34px;
  padding: 3px;
  display: flex;
}

.msp-canvas-mode-switch-segment {
  align-items: center;
  gap: var(--msp-space-2);
  height: 100%;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  transition: background-color var(--msp-motion-fast) var(--msp-ease),
    color var(--msp-motion-fast) var(--msp-ease);
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }

  &[aria-pressed="true"] {
    background: var(--msp-color-control);
    color: var(--msp-color-text);
  }
}

.msp-canvas-mode-switch-switch[data-canvas-mode="select"] .msp-canvas-mode-switch-segment[aria-pressed="true"] {
  background: color-mix(in oklab, var(--msp-color-accent), transparent 82%);
  color: var(--msp-color-accent);
}

.msp-canvas-mode-switch-level {
  animation: msp-canvas-mode-switch-level-in var(--msp-motion-base) var(--msp-ease);
  flex-direction: column;
  align-items: flex-end;
  display: flex;
  position: relative;
}

@keyframes msp-canvas-mode-switch-level-in {
  from {
    opacity: 0;
    transform: translateY(-4px);
  }
}

.msp-canvas-mode-switch-levelButton {
  align-items: center;
  gap: var(--msp-space-2);
  height: 30px;
  padding: 0 var(--msp-space-3);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  display: flex;

  &:hover, &[aria-expanded="true"] {
    background: var(--msp-glass-panel-bg);
  }
}

.msp-canvas-mode-switch-levelCaption {
  color: var(--msp-color-text-subtle);
}

.msp-canvas-mode-switch-chevron {
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="true"] > & {
    transform: rotate(180deg);
  }
}

.msp-canvas-mode-switch-levels {
  top: calc(100% + var(--msp-space-2));
  z-index: var(--msp-z-overlay);
  min-width: 168px;
  padding: var(--msp-space-2);
  background: var(--msp-glass-panel-bg);
  animation: msp-canvas-mode-switch-level-in var(--msp-motion-fast) var(--msp-ease);
  display: grid;
  position: absolute;
  right: 0;
}

.msp-canvas-mode-switch-group {
  padding: var(--msp-space-3) var(--msp-space-3) var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-canvas-mode-switch-levelRow {
  align-items: center;
  gap: var(--msp-space-3);
  height: 28px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  text-align: left;
  grid-template-columns: auto 1fr;
  display: grid;

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[aria-checked="true"] {
    color: var(--msp-color-accent);
    font-weight: 600;
  }
}

:where([data-mode="compact"]) .msp-canvas-mode-switch-anchor {
  right: var(--msp-canvas-rail-inset, var(--msp-space-4));
}

:where([data-mode="compact"]) .msp-canvas-mode-switch-label, :where([data-mode="compact"]) .msp-canvas-mode-switch-levelCaption {
  display: none;
}
}
`);
