/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from mode-chooser.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'card': 'msp-shell-mode-chooser-card',
  'description': 'msp-shell-mode-chooser-description',
  'heading': 'msp-shell-mode-chooser-heading',
  'icon': 'msp-shell-mode-chooser-icon',
  'note': 'msp-shell-mode-chooser-note',
  'option': 'msp-shell-mode-chooser-option',
  'options': 'msp-shell-mode-chooser-options',
  'rows': 'msp-shell-mode-chooser-rows',
  'title': 'msp-shell-mode-chooser-title',
  'welcome': 'msp-shell-mode-chooser-welcome',
} as const;

defineStyles('views/shell/mode-chooser.module.css', `@layer msp.components {
.msp-shell-mode-chooser-welcome {
  z-index: var(--msp-z-overlay);
  padding: var(--msp-space-5);
  pointer-events: auto;
  place-items: center;
  display: grid;
  position: absolute;
  inset: 0;
}

.msp-shell-mode-chooser-card {
  gap: var(--msp-space-4);
  width: min(560px, 100%);
  padding: var(--msp-space-5);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-panel-bg);
  backdrop-filter: blur(16px);
  box-shadow: 0 0 0 1px var(--msp-glass-edge),
    var(--msp-shadow-overlay);
  color: var(--msp-color-text);
  display: grid;
}

.msp-shell-mode-chooser-heading {
  font-size: var(--msp-font-size-lg);
  margin: 0;
  font-weight: 600;
}

.msp-shell-mode-chooser-note {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  margin: 0;
}

.msp-shell-mode-chooser-options {
  gap: var(--msp-space-3);
  grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
  display: grid;
}

.msp-shell-mode-chooser-option {
  align-items: center;
  column-gap: var(--msp-space-3);
  row-gap: var(--msp-space-1);
  padding: var(--msp-space-4);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  text-align: left;
  transition: background-color var(--msp-motion-fast) var(--msp-ease);
  grid-template-columns: auto minmax(0, 1fr);
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-checked="true"] {
    color: var(--msp-color-accent);
  }
}

.msp-shell-mode-chooser-icon {
  grid-row: span 2;
}

.msp-shell-mode-chooser-title {
  font-weight: 600;
  font-size: var(--msp-font-size-md);
}

.msp-shell-mode-chooser-description {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  grid-column: 2;
}

.msp-shell-mode-chooser-rows .msp-shell-mode-chooser-options {
  gap: var(--msp-space-1);
  grid-template-columns: minmax(0, 1fr);
}

.msp-shell-mode-chooser-rows .msp-shell-mode-chooser-option {
  padding: var(--msp-space-2) var(--msp-inset);
  background: none;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-shell-mode-chooser-rows .msp-shell-mode-chooser-title {
  font-size: var(--msp-font-size-sm);
}

.msp-shell-mode-chooser-rows .msp-shell-mode-chooser-description {
  font-size: var(--msp-font-size-xs);
}
}
`);
