/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from theme.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../registry.js';

defineStyles('styles/theme.css', `@layer msp.theme {
  :root, [data-msp-theme="light"] {
    color-scheme: light;
    --msp-color-bg: #eeece7;
    --msp-color-text: #332b1f;
    --msp-color-accent: #ae5d04;
    --msp-color-accent-strong: #974102;
    --msp-color-danger: #c4091f;
    --msp-color-success: #078a40;
    --msp-color-warning: #b7860b;
    --msp-color-info: #5e3673;
    --msp-shade: #000;
    --msp-color-brand-1: #7a4f2f;
    --msp-color-brand-2: #e0ae2c;
    --msp-color-brand-3: #6f9a80;
    --msp-color-brand-4: #3f86ad;
  }

  [data-msp-theme="dark"] {
    color-scheme: dark;
    --msp-color-bg: #111318;
    --msp-color-text: #ccd4e0;
    --msp-color-accent: #51a2fb;
    --msp-color-accent-strong: #68befd;
    --msp-color-danger: #fd354b;
    --msp-color-success: #0cca5d;
    --msp-color-warning: #fcc937;
    --msp-color-info: #b18ac6;
    --msp-shade: #fff;
    --msp-color-brand-1: #d29a68;
    --msp-color-brand-2: #f5d36b;
    --msp-color-brand-3: #a8d8bc;
    --msp-color-brand-4: #79bfe6;
  }

  [data-msp-theme="blue"] {
    color-scheme: dark;
    --msp-color-bg: #2d3e50;
    --msp-color-text: #edf1f2;
    --msp-color-accent: #3b9ad9;
    --msp-color-accent-strong: #fff;
    --msp-color-danger: #ff5a4f;
    --msp-color-success: #0cca5d;
    --msp-color-warning: #fcc937;
    --msp-color-info: #b18ac6;
    --msp-shade: #fff;
    --msp-color-brand-1: #d29a68;
    --msp-color-brand-2: #f5d36b;
    --msp-color-brand-3: #a8d8bc;
    --msp-color-brand-4: #79bfe6;
  }

  @media (prefers-color-scheme: dark) {
    [data-msp-theme="auto"] {
      color-scheme: dark;
      --msp-color-bg: #111318;
      --msp-color-text: #ccd4e0;
      --msp-color-accent: #51a2fb;
      --msp-color-accent-strong: #68befd;
      --msp-color-danger: #fd354b;
      --msp-color-success: #0cca5d;
      --msp-color-warning: #fcc937;
      --msp-color-info: #b18ac6;
      --msp-shade: #fff;
      --msp-color-brand-1: #d29a68;
      --msp-color-brand-2: #f5d36b;
      --msp-color-brand-3: #a8d8bc;
      --msp-color-brand-4: #79bfe6;
    }
  }

  :root, [data-msp-theme] {
    --msp-color-panel: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 3%);
    --msp-color-raised: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 5%);
    --msp-color-well: var(--msp-color-bg);
    --msp-color-control: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 6.5%);
    --msp-color-hover: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 9%);
    --msp-color-pressed: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 13%);
    --msp-color-border: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 15%);
    --msp-color-border-strong: color-mix(in oklab, var(--msp-color-bg), var(--msp-shade) 28%);
    --msp-color-text-muted: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 33%);
    --msp-color-text-subtle: color-mix(in oklab, var(--msp-color-text), var(--msp-color-bg) 55%);
    --msp-color-accent-muted: color-mix(in oklab, var(--msp-color-accent), var(--msp-color-bg) 80%);
    --msp-color-focus-ring: var(--msp-color-accent);
    --msp-color-scrim: #00000059;
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
    --msp-motion-base: .18s;
    --msp-ease: cubic-bezier(.2, 0, 0, 1);
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
