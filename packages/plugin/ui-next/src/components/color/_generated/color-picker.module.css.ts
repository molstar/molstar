/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from color-picker.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'custom': 'msp-color-color-picker-custom',
  'hex': 'msp-color-color-picker-hex',
  'palette': 'msp-color-color-picker-palette',
  'picker': 'msp-color-color-picker-picker',
  'preview': 'msp-color-color-picker-preview',
  'swatch': 'msp-color-color-picker-swatch',
} as const;

defineStyles('components/color/color-picker.module.css', `@layer msp.components {
.msp-color-color-picker-picker {
  gap: var(--msp-space-3);
  padding: var(--msp-space-3) var(--msp-space-3) var(--msp-space-3) var(--msp-space-4);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  box-shadow: inset 2px 0 0 var(--msp-color-accent);
  display: grid;
}

.msp-color-color-picker-palette {
  gap: var(--msp-space-2);
  grid-template-columns: repeat(8, minmax(0, 22px));
  justify-content: start;
  display: grid;
}

.msp-color-color-picker-swatch {
  aspect-ratio: 1;
  border-radius: var(--msp-radius-sm);
  min-height: 0;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  &:hover {
    transform: scale(1.08);
  }

  &[aria-selected="true"] {
    box-shadow: 0 0 0 1.5px var(--msp-color-bg),
      0 0 0 3px var(--msp-color-text-muted);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 2px;
  }
}

.msp-color-color-picker-custom {
  align-items: center;
  gap: var(--msp-space-2);
  display: flex;
}

.msp-color-color-picker-preview {
  border-radius: var(--msp-radius-sm);
  width: 22px;
  height: 22px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
  flex: none;
}

.msp-color-color-picker-hex {
  min-width: 0;
  height: 24px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  font-variant-numeric: tabular-nums;
  border: 0;
  flex: 1;

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 1px;
  }
}
}
`);
