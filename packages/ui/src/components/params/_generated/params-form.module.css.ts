/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from params-form.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'channels': 'msui-params-params-form-channels',
  'chevron': 'msui-params-params-form-chevron',
  'control': 'msui-params-params-form-control',
  'count': 'msui-params-params-form-count',
  'expand': 'msui-params-params-form-expand',
  'field': 'msui-params-params-form-field',
  'fields': 'msui-params-params-form-fields',
  'file': 'msui-params-params-form-file',
  'fileClear': 'msui-params-params-form-fileClear',
  'fileNames': 'msui-params-params-form-fileNames',
  'filePick': 'msui-params-params-form-filePick',
  'fold': 'msui-params-params-form-fold',
  'foldHeader': 'msui-params-params-form-foldHeader',
  'form': 'msui-params-params-form-form',
  'input': 'msui-params-params-form-input',
  'item': 'msui-params-params-form-item',
  'itemAction': 'msui-params-params-form-itemAction',
  'itemLabel': 'msui-params-params-form-itemLabel',
  'itemRow': 'msui-params-params-form-itemRow',
  'knob': 'msui-params-params-form-knob',
  'label': 'msui-params-params-form-label',
  'labelText': 'msui-params-params-form-labelText',
  'meta': 'msui-params-params-form-meta',
  'number': 'msui-params-params-form-number',
  'numeric': 'msui-params-params-form-numeric',
  'presets': 'msui-params-params-form-presets',
  'range': 'msui-params-params-form-range',
  'row': 'msui-params-params-form-row',
  'slider': 'msui-params-params-form-slider',
  'swatch': 'msui-params-params-form-swatch',
  'switch': 'msui-params-params-form-switch',
  'textarea': 'msui-params-params-form-textarea',
  'unsupported': 'msui-params-params-form-unsupported',
} as const;

defineStyles('components/params/params-form.module.css', `@layer msui.components {
.msui-params-params-form-form {
  --msui-params-level: 0;
  column-gap: var(--msui-space-3);
  padding-inline: var(--msui-inner);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  grid-template-columns: fit-content(50%) minmax(0, 1fr);
  display: grid;
}

.msui-params-params-form-field, .msui-params-params-form-fields, .msui-params-params-form-fold, .msui-params-params-form-item, .msui-params-params-form-row {
  grid-template-columns: subgrid;
  grid-column: 1 / -1;
  min-width: 0;
  display: grid;

  &[hidden] {
    display: none;
  }
}

:is(.msui-params-params-form-form, .msui-params-params-form-field, .msui-params-params-form-fields, .msui-params-params-form-fold, .msui-params-params-form-item) > * {
  grid-column: 1 / -1;
}

.msui-params-params-form-row {
  min-height: var(--msui-row);
  align-items: center;
}

.msui-params-params-form-row > .msui-params-params-form-label {
  grid-column: 1;
}

.msui-params-params-form-row > .msui-params-params-form-control {
  grid-column: 2;
}

.msui-params-params-form-label {
  min-width: var(--msui-label);
  padding-block: var(--msui-space-1);
  color: var(--msui-text-muted);
  white-space: nowrap;
  align-items: center;
  padding-inline-start: calc(var(--msui-params-level, 0) * var(--msui-form-indent));
  display: flex;
}

.msui-params-params-form-labelText {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-params-params-form-control {
  align-items: center;
  gap: var(--msui-space-1);
  min-width: 0;
  min-height: var(--msui-row);
  display: flex;
}

.msui-params-params-form-input, .msui-params-params-form-number, .msui-params-params-form-textarea {
  width: 100%;
  min-width: 0;
  height: var(--msui-control-height);
  padding: 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  color: var(--msui-text);
  font: inherit;
  border: 0;
  margin: 0;

  &:hover {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &::placeholder {
    color: var(--msui-text-subtle);
  }
}

.msui-params-params-form-number, .msui-params-params-form-numeric {
  font-variant-numeric: tabular-nums;
}

.msui-params-params-form-number {
  text-align: end;
}

.msui-params-params-form-textarea {
  height: auto;
  min-height: calc(3 * var(--msui-leading-sm) + 2 * var(--msui-space-1));
  padding-block: var(--msui-space-1);
  resize: vertical;
  display: block;
}

.msui-params-params-form-slider {
  align-items: center;
  gap: var(--msui-space-2);
  grid-template-columns: minmax(0, 1fr) 48px;
  width: 100%;
  display: grid;
}

.msui-params-params-form-range {
  --msui-params-fill: color-mix(in oklab, var(--msui-text-muted), transparent 45%);
  width: calc(100% - var(--msui-control-inset));
  height: var(--msui-control-height);
  margin: 0 0 0 var(--msui-control-inset);
  -webkit-appearance: none;
  appearance: none;
  cursor: pointer;
  background: none;

  &::-webkit-slider-runnable-track {
    border-radius: var(--msui-radius-sm);
    background: linear-gradient(to right,
      var(--msui-params-fill) var(--msui-params-progress, 0%),
      var(--msui-control-active) var(--msui-params-progress, 0%));
    height: 4px;
  }

  &::-webkit-slider-thumb {
    background: var(--msui-text-muted);
    -webkit-appearance: none;
    appearance: none;
    width: 12px;
    height: 12px;
    transition: transform var(--msui-motion-fast) var(--msui-ease);
    border-radius: 50%;
    margin-top: -4px;
  }

  &:hover::-webkit-slider-thumb, &:active::-webkit-slider-thumb {
    transform: scale(1.2);
  }

  &::-moz-range-track {
    border-radius: var(--msui-radius-sm);
    background: var(--msui-control-active);
    height: 4px;
  }

  &::-moz-range-progress {
    border-radius: var(--msui-radius-sm);
    background: var(--msui-params-fill);
    height: 4px;
  }

  &::-moz-range-thumb {
    background: var(--msui-text-muted);
    border: 0;
    border-radius: 50%;
    width: 12px;
    height: 12px;
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: 2px;
    border-radius: var(--msui-radius-sm);
  }
}

.msui-params-params-form-switch {
  width: 28px;
  height: 16px;
  background: var(--msui-control-active);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  border-radius: 8px;
  flex: none;
  margin-inline-start: var(--msui-control-inset);
  position: relative;

  &[aria-checked="true"] {
    background: var(--msui-text-muted);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: 2px;
  }
}

.msui-params-params-form-knob {
  background: var(--msui-bg);
  width: 12px;
  height: 12px;
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  border-radius: 50%;
  position: absolute;
  top: 2px;
  left: 2px;
  box-shadow: 0 1px 1px #0000001f;

  [aria-checked="true"] > & {
    transform: translateX(12px);
  }
}

.msui-params-params-form-swatch {
  border-radius: var(--msui-radius-sm);
  width: 12px;
  height: 12px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msui-text), transparent 88%);
  flex: none;
}

.msui-params-params-form-channels, .msui-params-params-form-fileNames, .msui-params-params-form-meta, .msui-params-params-form-itemLabel {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-params-params-form-channels {
  font-variant-numeric: tabular-nums;
}

.msui-params-params-form-unsupported {
  color: var(--msui-text-subtle);
  padding-inline-start: var(--msui-control-inset);
}

.msui-params-params-form-foldHeader {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: center;
  gap: var(--msui-space-1);
  min-width: 0;
  max-width: 100%;
  color: var(--msui-text-muted);
  text-align: start;
  display: flex;

  &:hover {
    color: var(--msui-text);
  }
}

.msui-params-params-form-chevron {
  color: var(--msui-text-subtle);
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  flex: none;

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}

.msui-params-params-form-meta {
  color: var(--msui-text-subtle);
  flex: 1;
  padding-inline-start: var(--msui-control-inset);
}

.msui-params-params-form-presets {
  flex: none;
  margin-inline-start: auto;
}

.msui-params-params-form-expand {
  flex: none;

  &[aria-expanded="true"] {
    background: var(--msui-control-active);
  }
}

.msui-params-params-form-switch ~ .msui-params-params-form-expand {
  margin-inline-start: auto;
}

.msui-params-params-form-file {
  width: 100%;
  min-width: 0;
  height: var(--msui-control-height);
  border-radius: var(--msui-radius-md);
  background: var(--msui-control);
  box-shadow: inset 0 0 0 1px var(--msui-field-edge);
  align-items: center;
  display: flex;

  &:hover {
    background: var(--msui-control-hover);
  }

  &[data-dragging] {
    background: var(--msui-control-hover);
    outline: 1px dashed var(--msui-line-strong);
    outline-offset: -1px;
  }
}

.msui-params-params-form-filePick {
  min-width: 0;
  height: 100%;
  padding: 0 var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  text-align: start;
  flex: 1;
  align-items: center;
  display: flex;

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-params-params-form-fileNames {
  color: var(--msui-text);
  flex: 1;

  &[data-empty] {
    color: var(--msui-text-muted);
  }
}

.msui-params-params-form-fileClear {
  --msui-icon-size: var(--msui-icon-sm);
  width: var(--msui-control-height);
  height: var(--msui-control-height);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  flex: none;
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msui-control-active);
    color: var(--msui-text);
  }
}

.msui-params-params-form-count {
  min-width: 0;
  color: var(--msui-text-subtle);
  flex: 1;
  padding-inline-start: var(--msui-control-inset);
}

.msui-params-params-form-item {
  display: grid;
}

.msui-params-params-form-itemRow {
  align-items: center;
  gap: var(--msui-space-1);
  min-height: var(--msui-row);
  margin-inline: calc(-1 * var(--msui-control-inset));
  padding-inline: calc(var(--msui-control-inset) + var(--msui-params-level, 0) * var(--msui-form-indent)) 0;
  border-radius: var(--msui-radius-md);
  cursor: pointer;
  display: flex;

  &:hover {
    background: var(--msui-control-hover);
  }
}

.msui-params-params-form-itemLabel {
  flex: 1;
}

.msui-params-params-form-itemRow > :last-child {
  margin-inline-end: var(--msui-control-inset);
}

.msui-params-params-form-itemAction {
  opacity: 0;
  flex: none;

  &[aria-expanded="true"] {
    background: var(--msui-control-active);
  }

  .msui-params-params-form-itemRow:hover > &, &:focus-visible, [data-open] > .msui-params-params-form-itemRow > & {
    opacity: 1;
  }
}

.msui-params-params-form-row:has( > .msui-params-params-form-control > .msui-params-params-form-textarea) {
  align-items: start;
}
}
`);
