/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from foundations.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../registry.js';

defineStyles('styles/foundations.css', `@layer msui.foundations {
  :root, :host, [data-msui-theme="light"] {
    color-scheme: light;
    --msui-bg: #f1efea;
    --msui-fg: #2d271e;
    --msui-accent: #8f4c06;
    --msui-danger: #b8101f;
    --msui-success: #12703a;
    --msui-warning: #7a5806;
    --msui-info: #5b3a72;
    --msui-on-accent: #fff;
    --msui-shade: #000;
    --msui-step-1: 3%;
    --msui-step-2: 6%;
    --msui-step-3: 9%;
    --msui-step-4: 13%;
    --msui-mix-muted: 25%;
    --msui-mix-subtle: 32%;
  }

  [data-msui-theme="dark"] {
    color-scheme: dark;
    --msui-bg: #121419;
    --msui-fg: #cfd6e1;
    --msui-accent: #5aa6f8;
    --msui-danger: #fb4a5c;
    --msui-success: #25c768;
    --msui-warning: #f2c43c;
    --msui-info: #b590c9;
    --msui-on-accent: #121419;
    --msui-shade: #fff;
    --msui-step-1: 4%;
    --msui-step-2: 9%;
    --msui-step-3: 14%;
    --msui-step-4: 19%;
    --msui-mix-muted: 26%;
    --msui-mix-subtle: 36%;
  }

  [data-msui-theme="blue"] {
    color-scheme: dark;
    --msui-bg: #2c3d4f;
    --msui-fg: #edf1f3;
    --msui-accent: #7cb9e6;
    --msui-danger: #ff9a8e;
    --msui-success: #3fd987;
    --msui-warning: #f6cf52;
    --msui-info: #c8aad6;
    --msui-on-accent: #15212d;
    --msui-shade: #fff;
    --msui-step-1: 4%;
    --msui-step-2: 9%;
    --msui-step-3: 14%;
    --msui-step-4: 19%;
    --msui-mix-muted: 18%;
    --msui-mix-subtle: 30%;
  }

  [data-msui-theme="contrast"] {
    color-scheme: dark;
    --msui-bg: #000;
    --msui-fg: #fff;
    --msui-accent: #ffd60a;
    --msui-danger: #ff6b6b;
    --msui-success: #3ddc84;
    --msui-warning: #ffd60a;
    --msui-info: #d0a8ff;
    --msui-on-accent: #000;
    --msui-shade: #fff;
    --msui-step-1: 8%;
    --msui-step-2: 16%;
    --msui-step-3: 24%;
    --msui-step-4: 32%;
    --msui-mix-muted: 12%;
    --msui-mix-subtle: 22%;
  }

  @media (prefers-color-scheme: dark) {
    [data-msui-theme="auto"] {
      color-scheme: dark;
      --msui-bg: #121419;
      --msui-fg: #cfd6e1;
      --msui-accent: #5aa6f8;
      --msui-danger: #fb4a5c;
      --msui-success: #25c768;
      --msui-warning: #f2c43c;
      --msui-info: #b590c9;
      --msui-on-accent: #121419;
      --msui-shade: #fff;
      --msui-step-1: 4%;
      --msui-step-2: 9%;
      --msui-step-3: 14%;
      --msui-step-4: 19%;
      --msui-mix-muted: 26%;
      --msui-mix-subtle: 36%;
    }
  }

  :root, :host, [data-msui-theme] {
    --msui-surface-sunken: var(--msui-bg);
    --msui-surface: color-mix(in oklab, var(--msui-bg), var(--msui-shade) var(--msui-step-1));
    --msui-surface-raised: color-mix(in oklab, var(--msui-bg), var(--msui-shade) var(--msui-step-2));
    --msui-control: var(--msui-surface-raised);
    --msui-control-hover: color-mix(in oklab, var(--msui-bg), var(--msui-shade) var(--msui-step-3));
    --msui-control-active: color-mix(in oklab, var(--msui-bg), var(--msui-shade) var(--msui-step-4));
    --msui-text: var(--msui-fg);
    --msui-text-muted: color-mix(in oklab, var(--msui-fg), var(--msui-bg) var(--msui-mix-muted));
    --msui-text-subtle: color-mix(in oklab, var(--msui-fg), var(--msui-bg) var(--msui-mix-subtle));
    --msui-line: color-mix(in oklab, var(--msui-bg), var(--msui-shade) 14%);
    --msui-line-strong: color-mix(in oklab, var(--msui-bg), var(--msui-shade) 28%);
    --msui-accent-muted: color-mix(in oklab, var(--msui-accent), var(--msui-bg) 80%);
    --msui-accent-strong: color-mix(in oklab, var(--msui-accent), var(--msui-shade) 12%);
    --msui-focus: var(--msui-accent);
    --msui-scrim: #00000059;
    --msui-field-edge: transparent;
    --msui-selected-edge: transparent;
    --msui-font-sans: "Helvetica Neue", "Segoe UI", Helvetica, "Source Sans Pro", Arial, sans-serif;
    --msui-font-mono: ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
    --msui-text-xs: 11px;
    --msui-leading-xs: 16px;
    --msui-text-sm: 12px;
    --msui-leading-sm: 16px;
    --msui-text-md: 14px;
    --msui-leading-md: 20px;
    --msui-text-lg: 16px;
    --msui-leading-lg: 24px;
    --msui-weight-regular: 400;
    --msui-weight-strong: 600;
    --msui-unit: 4px;
    --msui-space-1: calc(1 * var(--msui-unit));
    --msui-space-2: calc(2 * var(--msui-unit));
    --msui-space-3: calc(3 * var(--msui-unit));
    --msui-space-4: calc(4 * var(--msui-unit));
    --msui-space-5: calc(6 * var(--msui-unit));
    --msui-space-6: calc(8 * var(--msui-unit));
    --msui-row: calc(6 * var(--msui-unit));
    --msui-control-height: calc(var(--msui-row) - var(--msui-unit));
    --msui-header: calc(8 * var(--msui-unit));
    --msui-icon-sm: calc(3 * var(--msui-unit));
    --msui-icon: calc(4 * var(--msui-unit));
    --msui-icon-lg: calc(5 * var(--msui-unit));
    --msui-inset: var(--msui-space-3);
    --msui-surface-pad: var(--msui-space-1);
    --msui-inner: calc(var(--msui-inset) - var(--msui-surface-pad));
    --msui-icon-slot: var(--msui-icon);
    --msui-icon-gap: var(--msui-space-1);
    --msui-indent: calc(var(--msui-icon-slot) + var(--msui-icon-gap));
    --msui-form-indent: var(--msui-space-3);
    --msui-control-inset: var(--msui-space-2);
    --msui-label: calc(28 * var(--msui-unit));
    --msui-corner-shape: round;
    --msui-radius-sm: 2px;
    --msui-radius-lg: 8px;
    --msui-radius-md: max(0px, calc(var(--msui-radius-lg) - var(--msui-surface-pad)));
    --msui-edge: 0 0 0 1px var(--msui-line);
    --msui-surface-edge: none;
    --msui-shadow-float: 0 2px 8px #00000029;
    --msui-shadow-overlay: 0 8px 28px #00000042;
    --msui-z-overlay: 10;
    --msui-z-dialog: 20;
    --msui-z-toast: 30;
    --msui-z-tooltip: 40;
    --msui-motion-fast: .1s;
    --msui-motion-base: .2s;
    --msui-motion-slow: .36s;
    --msui-ease: cubic-bezier(.2, 0, 0, 1);
    --msui-ease-in-out: cubic-bezier(.4, 0, .2, 1);
    --msui-ease-in: cubic-bezier(.4, 0, 1, 1);
    --msui-enter: var(--msui-motion-base) var(--msui-ease);
    --msui-exit: var(--msui-motion-fast) var(--msui-ease-in);
    --msui-rise: var(--msui-space-1);
  }

  @supports (corner-shape: squircle) {
    :root, :host, [data-msui-theme] {
      --msui-corner-shape: squircle;
      --msui-radius-sm: 3px;
      --msui-radius-lg: 12px;
    }
  }

  [data-msui-theme="contrast"] {
    --msui-line: color-mix(in oklab, var(--msui-bg), var(--msui-shade) 45%);
    --msui-line-strong: color-mix(in oklab, var(--msui-bg), var(--msui-shade) 70%);
    --msui-surface-edge: inset 0 0 0 1px var(--msui-line-strong);
    --msui-field-edge: var(--msui-line-strong);
    --msui-selected-edge: var(--msui-accent);
  }

  @media (prefers-contrast: more) {
    [data-msui-theme] {
      --msui-text-muted: color-mix(in oklab, var(--msui-fg), var(--msui-bg) 15%);
      --msui-text-subtle: color-mix(in oklab, var(--msui-fg), var(--msui-bg) 25%);
      --msui-line: color-mix(in oklab, var(--msui-bg), var(--msui-shade) 40%);
      --msui-line-strong: color-mix(in oklab, var(--msui-bg), var(--msui-shade) 60%);
      --msui-surface-edge: inset 0 0 0 1px var(--msui-line-strong);
      --msui-field-edge: var(--msui-line-strong);
      --msui-selected-edge: var(--msui-accent);
    }
  }

  @media (prefers-reduced-motion: reduce) {
    :root, :host, [data-msui-theme] {
      --msui-motion-fast: 0s;
      --msui-motion-base: 0s;
      --msui-motion-slow: 0s;
    }
  }
}

@layer msui.primitives {
  .msui-surface {
    padding: var(--msui-space-1) var(--msui-surface-pad);
    border-radius: var(--msui-radius-lg);
    background: var(--msui-surface);
    box-shadow: var(--msui-surface-edge);
    color: var(--msui-text);
    font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
    align-content: start;
    display: grid;
  }

  .msui-header {
    align-items: center;
    gap: var(--msui-icon-gap);
    --msui-icon-size: var(--msui-icon);
    height: var(--msui-header);
    padding-inline: var(--msui-inner);
    font-weight: var(--msui-weight-strong);
    display: flex;
  }

  .msui-row {
    grid-template-columns: var(--msui-icon-slot) minmax(0, 1fr) auto;
    column-gap: var(--msui-icon-gap);
    height: var(--msui-row);
    border-radius: var(--msui-radius-md);
    align-items: center;
    padding-inline-start: calc(var(--msui-inner) + var(--msui-level, 0) * var(--msui-indent));
    padding-inline-end: var(--msui-inner);
    display: grid;
  }

  .msui-row-plain {
    grid-template-columns: minmax(0, 1fr) auto;
  }

  .msui-slot {
    --msui-icon-size: var(--msui-icon);
    width: var(--msui-icon-slot);
    height: var(--msui-icon);
    color: var(--msui-text-muted);
    justify-content: center;
    align-items: center;
    display: inline-flex;
  }

  .msui-form {
    grid-template-columns: var(--msui-label) minmax(0, 1fr);
    column-gap: var(--msui-space-3);
    grid-auto-rows: minmax(var(--msui-row), auto);
    padding-inline: var(--msui-inner);
    align-items: center;
    display: grid;
  }

  .msui-label {
    color: var(--msui-text-muted);
  }

  .msui-value {
    padding-inline-start: var(--msui-control-inset);
  }

  .msui-control {
    height: var(--msui-control-height);
    padding-inline: var(--msui-control-inset);
    border-radius: var(--msui-radius-md);
    background: var(--msui-control);
    color: var(--msui-text);
    font: inherit;
    border: 0;
    align-items: center;
    display: flex;
  }

  .msui-num {
    font-variant-numeric: tabular-nums;
    text-align: end;
  }
}
`);
