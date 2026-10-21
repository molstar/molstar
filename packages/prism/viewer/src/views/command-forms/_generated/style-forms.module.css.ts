/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from style-forms.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'caption': 'msp-viewer-style-forms-caption',
  'choice': 'msp-viewer-style-forms-choice',
  'choiceLabel': 'msp-viewer-style-forms-choiceLabel',
  'empty': 'msp-viewer-style-forms-empty',
  'form': 'msp-viewer-style-forms-form',
  'grid': 'msp-viewer-style-forms-grid',
  'icon': 'msp-viewer-style-forms-icon',
  'label': 'msp-viewer-style-forms-label',
  'meta': 'msp-viewer-style-forms-meta',
  'one': 'msp-viewer-style-forms-one',
  'options': 'msp-viewer-style-forms-options',
  'swatch': 'msp-viewer-style-forms-swatch',
  'tile': 'msp-viewer-style-forms-tile',
  'unpainted': 'msp-viewer-style-forms-unpainted',
} as const;

defineStyles('viewer/views/command-forms/style-forms.module.css', `@layer msp.components {
.msp-viewer-style-forms-form {
  gap: var(--msp-space-1);
  padding: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-style-forms-caption {
  padding: var(--msp-space-1) var(--msp-space-3) 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-style-forms-grid {
  grid-template-columns: repeat(var(--columns, 3), minmax(0, 1fr));
  gap: var(--msp-space-1);
  display: grid;
}

@media (width <= 480px) {
  .msp-viewer-style-forms-grid {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}

.msp-viewer-style-forms-tile {
  align-items: center;
  gap: var(--msp-space-2);
  height: 28px;
  padding: 0 var(--msp-inset);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text);
  text-align: left;
  grid-template-columns: auto minmax(0, 1fr) auto;
  display: grid;

  &:hover, &[data-highlight="true"] {
    background: var(--msp-color-hover);
  }

  &[aria-checked="true"] {
    background: var(--msp-color-control);
    box-shadow: inset 0 0 0 1px var(--msp-selected-edge);
  }

  &[aria-checked="true"] .msp-viewer-style-forms-icon, &[aria-checked="true"] .msp-viewer-style-forms-label {
    color: var(--msp-color-accent);
  }
}

.msp-viewer-style-forms-icon {
  color: var(--msp-color-text-muted);
}

.msp-viewer-style-forms-label {
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-style-forms-meta {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-style-forms-swatch {
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  width: 24px;
  height: 14px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
  position: relative;
}

.msp-viewer-style-forms-unpainted {
  background: repeating-linear-gradient(135deg,
    var(--msp-color-control) 0 3px,
    color-mix(in oklab, var(--msp-color-text-muted), transparent 75%) 3px 6px);
}

.msp-viewer-style-forms-one {
  background: var(--msp-color-bg);
  width: 4px;
  height: 4px;
  box-shadow: 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 70%);
  border-radius: 50%;
  margin: -2px 0 0 -2px;
  position: absolute;
  top: 50%;
  left: 50%;
}

.msp-viewer-style-forms-empty {
  padding: var(--msp-space-3);
  color: var(--msp-color-text-subtle);
  text-align: center;
  grid-column: 1 / -1;
  margin: 0;
}

.msp-viewer-style-forms-options {
  gap: var(--msp-space-2);
  margin-top: var(--msp-space-2);
  padding: var(--msp-space-3) var(--msp-space-3) var(--msp-space-1);
  border-top: 1px solid var(--msp-glass-edge);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-style-forms-choice {
  align-items: center;
  gap: var(--msp-space-2);
  grid-template-columns: 72px minmax(0, 1fr);
  display: grid;
}

.msp-viewer-style-forms-choiceLabel {
  color: var(--msp-color-text-muted);
}
}
`);
