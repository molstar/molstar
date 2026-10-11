/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from mode-switch.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'granularity': 'msp-canvas-mode-switch-granularity',
  'granularity-in': 'msp-canvas-mode-switch-granularity-in',
  'label': 'msp-canvas-mode-switch-label',
  'segment': 'msp-canvas-mode-switch-segment',
  'switch': 'msp-canvas-mode-switch-switch',
} as const;

defineStyles('views/canvas/mode-switch.module.css', `@layer msp.components {
.msp-canvas-mode-switch-switch {
  top: var(--msp-canvas-top-inset, var(--msp-space-4));
  right: var(--msp-canvas-right-inset, var(--msp-space-4));
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  height: 34px;
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  transition: right var(--msp-motion-base) var(--msp-ease);
  align-items: center;
  gap: 2px;
  padding: 3px;
  display: flex;
  position: absolute;
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

.msp-canvas-mode-switch-granularity {
  height: 100%;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  animation: msp-canvas-mode-switch-granularity-in var(--msp-motion-base) var(--msp-ease);
  background: none;
  border: 0;

  &:hover {
    background: var(--msp-color-hover);
  }
}

@keyframes msp-canvas-mode-switch-granularity-in {
  from {
    opacity: 0;
    transform: translateX(-4px);
  }
}

:where([data-mode="compact"]) .msp-canvas-mode-switch-switch {
  right: var(--msp-canvas-rail-inset, var(--msp-space-4));
}

:where([data-mode="compact"]) .msp-canvas-mode-switch-label {
  display: none;
}
}
`);
