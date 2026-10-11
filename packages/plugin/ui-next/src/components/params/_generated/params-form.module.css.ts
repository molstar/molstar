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
  'range': 'msp-params-params-form-range',
  'row': 'msp-params-params-form-row',
  'select': 'msp-params-params-form-select',
  'settings': 'msp-params-params-form-settings',
  'slider': 'msp-params-params-form-slider',
  'subSettings': 'msp-params-params-form-subSettings',
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

.msp-params-params-form-input, .msp-params-params-form-number {
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
  width: 100%;
}

.msp-params-params-form-slider {
  align-items: center;
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr) 52px;
  width: 100%;
  display: grid;
}

.msp-params-params-form-range {
  -webkit-appearance: none;
  appearance: none;
  cursor: pointer;
  background: none;
  width: 100%;
  height: 16px;
  margin: 0;

  &::-webkit-slider-runnable-track {
    background: linear-gradient(to right,
      var(--msp-color-text-muted) var(--msp-progress, 0%),
      var(--msp-color-border-strong) var(--msp-progress, 0%));
    border-radius: 2px;
    height: 4px;
  }

  &::-webkit-slider-thumb {
    background: var(--msp-color-text);
    -webkit-appearance: none;
    appearance: none;
    width: 12px;
    height: 12px;
    transition: transform var(--msp-motion-fast) var(--msp-ease);
    border-radius: 50%;
    margin-top: -4px;
  }

  &:hover::-webkit-slider-thumb, &:active::-webkit-slider-thumb {
    transform: scale(1.25);
  }

  &::-moz-range-track {
    background: var(--msp-color-border-strong);
    border-radius: 2px;
    height: 4px;
  }

  &::-moz-range-progress {
    background: var(--msp-color-text-muted);
    border-radius: 2px;
    height: 4px;
  }

  &::-moz-range-thumb {
    background: var(--msp-color-text);
    border: 0;
    border-radius: 50%;
    width: 12px;
    height: 12px;
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 2px;
    border-radius: var(--msp-radius-sm);
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
    background: var(--msp-color-text-muted);
  }

  &:focus-visible {
    outline: 2px solid var(--msp-color-focus-ring);
    outline-offset: 2px;
  }
}

.msp-params-params-form-knob {
  background: var(--msp-color-bg);
  width: 14px;
  height: 14px;
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  border-radius: 50%;
  position: absolute;
  top: 2px;
  left: 2px;
  box-shadow: 0 1px 1px #0000001f;

  [aria-checked="true"] > & {
    transform: translateX(14px);
  }
}

.msp-params-params-form-color {
  align-items: center;
  gap: var(--msp-space-2);
  width: 100%;
  height: 24px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  display: flex;

  &:hover, &[aria-expanded="true"] {
    background: var(--msp-color-hover);
  }
}

.msp-params-params-form-swatch {
  border-radius: var(--msp-radius-sm);
  width: 16px;
  height: 16px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
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

.msp-params-params-form-settings {
  border-radius: var(--msp-radius-md);
  width: 24px;
  height: 24px;
  color: var(--msp-color-text-muted);
  place-items: center;
  margin-left: auto;
  display: grid;

  &:hover, &[aria-expanded="true"] {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-params-params-form-subSettings {
  margin: 2px 0 var(--msp-space-1);
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-1) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  box-shadow: inset 2px 0 0 var(--msp-color-accent);
}
}
`);
