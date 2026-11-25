/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from color-picker.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'custom': 'msui-color-color-picker-custom',
  'field': 'msui-color-color-picker-field',
  'greys': 'msui-color-color-picker-greys',
  'marker': 'msui-color-color-picker-marker',
  'palette': 'msui-color-color-picker-palette',
  'picker': 'msui-color-color-picker-picker',
  'preview': 'msui-color-color-picker-preview',
  'swatch': 'msui-color-color-picker-swatch',
  'text': 'msui-color-color-picker-text',
} as const;

defineStyles('components/color/color-picker.module.css', `@layer msui.components {
.msui-color-color-picker-picker {
  gap: var(--msui-space-2);
  display: grid;
}

.msui-color-color-picker-field, .msui-color-color-picker-greys {
  border-radius: var(--msui-radius-sm);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  cursor: crosshair;
  touch-action: none;
  position: relative;
}

.msui-color-color-picker-field {
  background: linear-gradient(#ffffffe6, #fff0 50%, #0000 50%, #000000e6), linear-gradient(to right, red, #ff0, #0f0, #0ff, #00f, #f0f, red);
  height: 96px;
}

.msui-color-color-picker-greys {
  background: linear-gradient(to right, #000, #fff);
  height: 12px;
}

.msui-color-color-picker-marker {
  pointer-events: none;
  border-radius: 50%;
  width: 10px;
  height: 10px;
  position: absolute;
  translate: -50% -50%;
  box-shadow: 0 0 0 1.5px #fff, 0 0 0 2.5px #00000080;
}

.msui-color-color-picker-palette {
  gap: var(--msui-space-1);
  grid-template-columns: repeat(8, minmax(0, 1fr));
  display: grid;
}

.msui-color-color-picker-swatch {
  border-radius: var(--msui-radius-sm);
  min-width: 0;
  height: 16px;
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  transition: transform var(--msui-motion-fast) var(--msui-ease);

  &[data-tone="edge"] {
    box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msui-text), transparent 85%);
  }

  &:hover {
    transform: scale(1.08);
  }

  &[aria-checked="true"] {
    outline: 2px solid var(--msui-text-muted);
    outline-offset: 1px;
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: 2px;
  }
}

.msui-color-color-picker-custom {
  align-items: center;
  gap: var(--msui-space-1);
  display: flex;
}

.msui-color-color-picker-preview {
  width: var(--msui-control-height);
  height: var(--msui-control-height);
  border-radius: var(--msui-radius-sm);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  flex: none;

  &[data-tone="edge"] {
    box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msui-text), transparent 85%);
  }
}

.msui-color-color-picker-text {
  min-width: 0;
  height: var(--msui-control-height);
  padding: 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  color: var(--msui-text);
  font: inherit;
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  font-variant-numeric: tabular-nums;
  border: 0;
  flex: 1;
  margin: 0;

  &::placeholder {
    color: var(--msui-text-subtle);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}
}
`);
