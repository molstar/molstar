/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from params-form.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'chevron': 'msp-params-params-form-chevron',
  'collapsible': 'msp-params-params-form-collapsible',
  'collapsibleHeader': 'msp-params-params-form-collapsibleHeader',
  'color': 'msp-params-params-form-color',
  'control': 'msp-params-params-form-control',
  'form': 'msp-params-params-form-form',
  'hex': 'msp-params-params-form-hex',
  'input': 'msp-params-params-form-input',
  'knob': 'msp-params-params-form-knob',
  'label': 'msp-params-params-form-label',
  'mapped': 'msp-params-params-form-mapped',
  'nested': 'msp-params-params-form-nested',
  'number': 'msp-params-params-form-number',
  'row': 'msp-params-params-form-row',
  'select': 'msp-params-params-form-select',
  'slider': 'msp-params-params-form-slider',
  'swatch': 'msp-params-params-form-swatch',
  'switch': 'msp-params-params-form-switch',
  'unsupported': 'msp-params-params-form-unsupported',
} as const;

defineStyles('components/params/params-form.module.css', `@layer msp.components {
.msp-params-params-form-form {
  gap: 2px;
  display: grid;
}

.msp-params-params-form-row {
  grid-template-columns: var(--msp-label-width) minmax(0, 1fr);
  align-items: center;
  gap: var(--msp-space-3);
  min-height: 24px;
  display: grid;

  &[hidden] {
    display: none;
  }
}

.msp-params-params-form-label {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  text-overflow: ellipsis;
  white-space: nowrap;
  overflow: hidden;
}

.msp-params-params-form-control {
  align-items: center;
  min-width: 0;
  display: flex;
}

.msp-params-params-form-input, .msp-params-params-form-select, .msp-params-params-form-number {
  width: 100%;
  height: 24px;
  padding: 0 var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  border: 0;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-params-params-form-select {
  padding-right: var(--msp-space-2);
}

.msp-params-params-form-slider {
  align-items: center;
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr) 52px;
  width: 100%;
  display: grid;

  & input[type="range"] {
    width: 100%;
    accent-color: var(--msp-color-accent);
  }
}

.msp-params-params-form-number {
  padding: 0 var(--msp-space-2);
  font-variant-numeric: tabular-nums;
  text-align: right;
}

.msp-params-params-form-switch {
  background: var(--msp-color-pressed);
  width: 32px;
  height: 18px;
  transition: background-color var(--msp-motion-fast) var(--msp-ease);
  border-radius: 9px;
  position: relative;

  &[aria-checked="true"] {
    background: var(--msp-color-accent);
  }
}

.msp-params-params-form-knob {
  width: 14px;
  height: 14px;
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  background: #fff;
  border-radius: 50%;
  position: absolute;
  top: 2px;
  left: 2px;

  [aria-checked="true"] > & {
    transform: translateX(14px);
  }
}

.msp-params-params-form-color {
  align-items: center;
  gap: var(--msp-space-2);
  height: 24px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  cursor: pointer;
  display: flex;

  & input[type="color"] {
    opacity: 0;
    width: 0;
    height: 0;
  }
}

.msp-params-params-form-swatch {
  border-radius: 50%;
  width: 16px;
  height: 16px;
}

.msp-params-params-form-hex {
  color: var(--msp-color-text-muted);
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
}

.msp-params-params-form-unsupported {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-params-params-form-mapped {
  gap: 2px;
  display: grid;
}

.msp-params-params-form-nested {
  padding-left: var(--msp-space-3);
  gap: 2px;
  display: grid;
}

.msp-params-params-form-collapsible {
  display: grid;
}

.msp-params-params-form-collapsibleHeader {
  align-items: center;
  gap: var(--msp-space-1);
  height: 28px;
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  font-weight: 600;
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-params-params-form-chevron {
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}
}
`);
