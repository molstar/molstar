/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from theme.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../registry.js';

defineStyles('styles/theme.css', `@layer msp.theme {
  :root, :host, [data-msp-theme="light"] {
    color-scheme: light;
    --msp-color-bg: #eeece7;
    --msp-canvas-bg: #fcfbf9;
    --msp-color-text: #332b1f;
    --msp-color-accent: #995103;
    --msp-color-accent-strong: #974102;
    --msp-color-danger: #c4091f;
    --msp-color-success: #057434;
    --msp-color-warning: #805d05;
    --msp-color-info: #5e3673;
    --msp-color-on-accent: #fff;
    --msp-shade: #000;
    --msp-mix-muted: 24%;
    --msp-mix-subtle: 31%;
    --msp-color-brand-1: #7a4f2f;
    --msp-color-brand-2: #e0ae2c;
    --msp-color-brand-3: #6f9a80;
    --msp-color-brand-4: #3f86ad;
  }

  [data-msp-theme="dark"] {
    color-scheme: dark;
    --msp-color-bg: #111318;
    --msp-canvas-bg: #16181d;
    --msp-color-text: #ccd4e0;
    --msp-color-accent: #51a2fb;
    --msp-color-accent-strong: #68befd;
    --msp-color-danger: #fd354b;
    --msp-color-success: #0cca5d;
    --msp-color-warning: #fcc937;
    --msp-color-info: #b18ac6;
    --msp-color-on-accent: #111318;
    --msp-shade: #fff;
    --msp-mix-muted: 25%;
    --msp-mix-subtle: 35%;
    --msp-step-control: 10%;
    --msp-step-hover: 15%;
    --msp-step-pressed: 20%;
    --msp-color-brand-1: #d29a68;
    --msp-color-brand-2: #f5d36b;
    --msp-color-brand-3: #a8d8bc;
    --msp-color-brand-4: #79bfe6;
  }

  [data-msp-theme="blue"] {
    color-scheme: dark;
    --msp-color-bg: #2d3e50;
    --msp-canvas-bg: #243344;
    --msp-color-text: #edf1f2;
    --msp-color-accent: #77b6e5;
    --msp-color-accent-strong: #fff;
    --msp-color-danger: #ff9083;
    --msp-color-success: #0cca5d;
    --msp-color-warning: #fcc937;
    --msp-color-info: #c3a4d3;
    --msp-color-on-accent: #14202c;
    --msp-shade: #fff;
    --msp-mix-muted: 18%;
    --msp-mix-subtle: 30%;
    --msp-step-control: 10%;
    --msp-step-hover: 15%;
    --msp-step-pressed: 20%;
    --msp-color-brand-1: #d29a68;
    --msp-color-brand-2: #f5d36b;
    --msp-color-brand-3: #a8d8bc;
    --msp-color-brand-4: #79bfe6;
  }

  [data-msp-theme="contrast"] {
    color-scheme: dark;
    --msp-color-bg: #000;
    --msp-canvas-bg: #000;
    --msp-color-text: #fff;
    --msp-color-accent: #ffd60a;
    --msp-color-accent-strong: #ffe866;
    --msp-color-danger: #ff6b6b;
    --msp-color-success: #3ddc84;
    --msp-color-warning: #ffd60a;
    --msp-color-info: #d0a8ff;
    --msp-color-on-accent: #000;
    --msp-shade: #fff;
    --msp-step-control: 16%;
    --msp-step-hover: 24%;
    --msp-step-pressed: 32%;
    --msp-color-brand-1: #ffb37a;
    --msp-color-brand-2: #ffe066;
    --msp-color-brand-3: #9ff0c4;
    --msp-color-brand-4: #8fd3ff;
  }

  @media (prefers-color-scheme: dark) {
    [data-msp-theme="auto"] {
      color-scheme: dark;
      --msp-color-bg: #111318;
      --msp-canvas-bg: #16181d;
      --msp-color-text: #ccd4e0;
      --msp-color-accent: #51a2fb;
      --msp-color-accent-strong: #68befd;
      --msp-color-danger: #fd354b;
      --msp-color-success: #0cca5d;
      --msp-color-warning: #fcc937;
      --msp-color-info: #b18ac6;
      --msp-color-on-accent: #111318;
      --msp-shade: #fff;
      --msp-mix-muted: 25%;
      --msp-mix-subtle: 35%;
      --msp-step-control: 10%;
      --msp-step-hover: 15%;
      --msp-step-pressed: 20%;
      --msp-color-brand-1: #d29a68;
      --msp-color-brand-2: #f5d36b;
      --msp-color-brand-3: #a8d8bc;
      --msp-color-brand-4: #79bfe6;
    }
  }

  :root, :host, [data-msp-theme] {
    --msp-color-panel: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 3%);
    --msp-color-raised: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 5%);
    --msp-color-well: var(--msp-color-bg);
    --msp-color-control: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) var(--msp-step-control, 6.5%));
    --msp-color-hover: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) var(--msp-step-hover, 9%));
    --msp-color-pressed: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) var(--msp-step-pressed, 13%));
    --msp-color-border: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 15%);
    --msp-color-border-strong: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 28%);
    --msp-color-text-muted: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) var(--msp-mix-muted, 24%));
    --msp-color-text-subtle: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) var(--msp-mix-subtle, 31%));
    --msp-color-accent-muted: color-mix(in oklab, var(--msp-color-accent), var(--msp-color-bg) 80%);
    --msp-color-focus-ring: var(--msp-color-accent);
    --msp-color-scrim: #00000059;
    --msp-field-edge: transparent;
    --msp-selected-edge: transparent;
    --msp-glass-bg: color-mix(in oklab, var(--msp-color-panel), transparent 30%);
    --msp-glass-panel-bg: color-mix(in oklab, var(--msp-color-panel), transparent 5%);
    --msp-glass-filter: blur(20px) saturate(1.4);
    --msp-glass-edge: color-mix(in oklab, var(--msp-color-text), transparent 88%);
    --msp-glass-shadow: 0 0 0 1px var(--msp-glass-edge);
    --msp-accent-cyan: #7f8c8d;
    --msp-accent-red: #e74c3c;
    --msp-accent-gray: #34495e;
    --msp-accent-green: #27ae60;
    --msp-accent-purple: #8e44ad;
    --msp-accent-blue: #2e86c1;
    --msp-accent-orange: #d35400;
    --msp-font-sans: "Helvetica Neue", "Segoe UI", Helvetica, "Source Sans Pro", Arial, sans-serif;
    --msp-font-mono: ui-monospace, "SF Mono", Menlo, Consolas, "Liberation Mono", monospace;
    --msp-font-size-xs: 11px;
    --msp-font-size-sm: 12px;
    --msp-font-size-md: 14px;
    --msp-font-size-lg: 16px;
    --msp-line-height: 1.4;
    --msp-row-height: 32px;
    --msp-control-height: 28px;
    --msp-label-width: 110px;
    --msp-space-1: 2px;
    --msp-space-2: 4px;
    --msp-space-3: 8px;
    --msp-space-4: 12px;
    --msp-space-5: 16px;
    --msp-space-6: 24px;
    --msp-inset: var(--msp-space-3);
    --msp-panel-inset: var(--msp-space-2);
    --msp-canvas-edge: var(--msp-space-3);
    --msp-icon-sm: 14px;
    --msp-icon-md: 15px;
    --msp-radius-sm: 2px;
    --msp-radius-md: 4px;
    --msp-radius-lg: 8px;
    --msp-shadow-float: 0 2px 8px #0000002e;
    --msp-shadow-overlay: 0 8px 28px #00000047;
    --msp-z-overlay: 10;
    --msp-z-dialog: 20;
    --msp-z-toast: 30;
    --msp-z-tooltip: 40;
    --msp-motion-fast: .1s;
    --msp-motion-base: .22s;
    --msp-motion-fade: .36s;
    --msp-motion-reveal: .9s;
    --msp-ease: cubic-bezier(.2, 0, 0, 1);
    --msp-ease-reveal: cubic-bezier(.4, 0, .2, 1);
  }

  [data-msp-theme="contrast"] {
    --msp-color-text-muted: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 12%);
    --msp-color-text-subtle: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 25%);
    --msp-color-border: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 45%);
    --msp-color-border-strong: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 70%);
    --msp-color-accent-muted: color-mix(in oklab, var(--msp-color-accent), var(--msp-color-bg) 60%);
    --msp-glass-bg: var(--msp-color-panel);
    --msp-glass-panel-bg: var(--msp-color-panel);
    --msp-glass-edge: var(--msp-color-border-strong);
    --msp-field-edge: var(--msp-color-border-strong);
    --msp-selected-edge: var(--msp-color-accent);
  }

  @media (prefers-contrast: more) {
    [data-msp-theme] {
      --msp-color-text-muted: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 15%);
      --msp-color-text-subtle: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 30%);
      --msp-color-border: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 40%);
      --msp-color-border-strong: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 60%);
      --msp-glass-bg: var(--msp-color-panel);
      --msp-glass-panel-bg: var(--msp-color-panel);
      --msp-glass-edge: var(--msp-color-border-strong);
      --msp-field-edge: var(--msp-color-border-strong);
      --msp-selected-edge: var(--msp-color-accent);
    }
  }

  [data-msp-ink="light"] {
    --msp-ink: #f5f5f5;
    --msp-ink-muted: #f5f5f5b3;
    --msp-ink-shadow: 0 0 3px #000c;
  }

  [data-msp-ink="dark"] {
    --msp-ink: #1a1a1a;
    --msp-ink-muted: #1a1a1ab3;
    --msp-ink-shadow: 0 0 3px #fffc;
  }
}
`);
