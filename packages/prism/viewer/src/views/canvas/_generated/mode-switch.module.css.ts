/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from mode-switch.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'badge': 'msp-viewer-mode-switch-badge',
  'chevron': 'msp-viewer-mode-switch-chevron',
  'group': 'msp-viewer-mode-switch-group',
  'level': 'msp-viewer-mode-switch-level',
  'level-in': 'msp-viewer-mode-switch-level-in',
  'levelButton': 'msp-viewer-mode-switch-levelButton',
  'levelCaption': 'msp-viewer-mode-switch-levelCaption',
  'levelRow': 'msp-viewer-mode-switch-levelRow',
  'levels': 'msp-viewer-mode-switch-levels',
  'segment': 'msp-viewer-mode-switch-segment',
  'select': 'msp-viewer-mode-switch-select',
  'switch': 'msp-viewer-mode-switch-switch',
} as const;

defineStyles('viewer/views/canvas/mode-switch.module.css', `@layer msp.components {
.msp-viewer-mode-switch-levelButton, .msp-viewer-mode-switch-levels {
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
}

.msp-viewer-mode-switch-switch {
  flex-direction: column;
  display: flex;
}

.msp-viewer-mode-switch-segment {
  width: 100%;
  height: var(--msp-stack-row, 40px);
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
    color: var(--msp-color-text);
  }
}

.msp-viewer-mode-switch-switch[data-canvas-mode="select"] .msp-viewer-mode-switch-segment[aria-pressed="true"] {
  background: color-mix(in oklab, var(--msp-color-accent), transparent 82%);
  color: var(--msp-color-accent);
}

.msp-viewer-mode-switch-level {
  top: var(--msp-stack-row, 40px);
  padding-right: var(--msp-space-3);
  opacity: 0;
  visibility: hidden;
  transition: opacity var(--msp-motion-fast) var(--msp-ease),
    transform var(--msp-motion-fast) var(--msp-ease),
    visibility 0s linear .16s;
  flex-direction: column;
  align-items: flex-end;
  display: flex;
  position: absolute;
  right: 100%;
  transform: translateX(4px);
}

:where(:has( > * .msp-viewer-mode-switch-select:hover,  > * .msp-viewer-mode-switch-select:focus-visible)) > .msp-viewer-mode-switch-level, .msp-viewer-mode-switch-level:hover, .msp-viewer-mode-switch-level:focus-within, .msp-viewer-mode-switch-level:has([aria-expanded="true"]) {
  opacity: 1;
  visibility: visible;
  transition-delay: 0s;
  transform: none;
}

@media (hover: none) {
  .msp-viewer-mode-switch-level {
    opacity: 1;
    visibility: visible;
    transform: none;
  }
}

.msp-viewer-mode-switch-select {
  position: relative;
}

.msp-viewer-mode-switch-badge {
  opacity: .8;
  position: absolute;
  bottom: 4px;
  right: 4px;
}

.msp-viewer-mode-switch-levelButton {
  align-items: center;
  gap: var(--msp-space-2);
  height: var(--msp-stack-row, 40px);
  padding: 0 var(--msp-inset);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  white-space: nowrap;
  display: flex;

  &:hover, &[aria-expanded="true"] {
    background: var(--msp-glass-panel-bg);
  }
}

.msp-viewer-mode-switch-levelCaption {
  color: var(--msp-color-text-subtle);
}

.msp-viewer-mode-switch-chevron {
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="true"] > & {
    transform: rotate(180deg);
  }
}

.msp-viewer-mode-switch-levels {
  top: calc(100% + var(--msp-space-2));
  right: var(--msp-space-3);
  z-index: var(--msp-z-overlay);
  min-width: 168px;
  padding: var(--msp-space-2);
  background: var(--msp-glass-panel-bg);
  animation: msp-viewer-mode-switch-level-in var(--msp-motion-fast) var(--msp-ease);
  display: grid;
  position: absolute;
}

.msp-viewer-mode-switch-group {
  padding: var(--msp-space-3) var(--msp-space-3) var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-mode-switch-levelRow {
  align-items: center;
  gap: var(--msp-space-3);
  height: 28px;
  padding: 0 var(--msp-inset);
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

:where([data-mode="compact"]) .msp-viewer-mode-switch-levelCaption {
  display: none;
}
}
`);
