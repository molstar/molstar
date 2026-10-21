/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from selection-bar.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'action': 'msp-viewer-selection-bar-action',
  'bar': 'msp-viewer-selection-bar-bar',
  'bar-in': 'msp-viewer-selection-bar-bar-in',
  'menu': 'msp-viewer-selection-bar-menu',
  'menuCard': 'msp-viewer-selection-bar-menuCard',
  'menuMark': 'msp-viewer-selection-bar-menuMark',
  'summary': 'msp-viewer-selection-bar-summary',
} as const;

defineStyles('viewer/views/canvas/selection-bar.module.css', `@layer msp.components {
.msp-viewer-selection-bar-bar {
  height: max(34px, var(--msp-touch-target, 0px));
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  animation: msp-viewer-selection-bar-bar-in var(--msp-motion-base) var(--msp-ease);
  align-items: center;
  padding: 0;
  display: flex;
  position: relative;
}

@keyframes msp-viewer-selection-bar-bar-in {
  from {
    opacity: 0;
    transform: translateY(4px);
  }
}

.msp-viewer-selection-bar-summary {
  align-items: center;
  gap: var(--msp-space-2);
  height: 100%;
  padding: 0 var(--msp-space-4);
  border-radius: var(--msp-radius-lg) 0 0 var(--msp-radius-lg);
  font-size: var(--msp-font-size-sm);
  white-space: nowrap;
  display: flex;
  position: relative;

  &:hover, &[aria-expanded="true"] {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-selection-bar-menuMark {
  opacity: .45;
  border-top: 4px solid #0000;
  border-right: 4px solid;
  position: absolute;
  bottom: 5px;
  right: 5px;
}

.msp-viewer-selection-bar-action {
  width: max(34px, var(--msp-touch-target, 0px));
  border-radius: 0 var(--msp-radius-lg) var(--msp-radius-lg) 0;
  height: 100%;
  color: var(--msp-color-text-muted);
  flex: none;
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-viewer-selection-bar-menu {
  z-index: var(--msp-z-overlay);
  padding-bottom: var(--msp-space-3);
  animation: msp-viewer-selection-bar-bar-in var(--msp-motion-fast) var(--msp-ease);
  position: absolute;
  bottom: 100%;
  right: 0;
}

.msp-viewer-selection-bar-menuCard {
  min-width: 200px;
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-panel-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  white-space: nowrap;
  display: grid;
}
}
`);
