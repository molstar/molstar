/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from glass-wordmark.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'canvas': 'msui-logo-glass-wordmark-canvas',
  'probe': 'msui-logo-glass-wordmark-probe',
  'wordmark': 'msui-logo-glass-wordmark-wordmark',
} as const;

defineStyles('components/logo/glass-wordmark.module.css', `@layer msui.components {
.msui-logo-glass-wordmark-wordmark {
  vertical-align: middle;
  -webkit-user-select: none;
  user-select: none;
  -webkit-tap-highlight-color: transparent;
  flex: none;
  display: inline-block;
  position: relative;
}

.msui-logo-glass-wordmark-canvas {
  pointer-events: none;
  position: absolute;
}

.msui-logo-glass-wordmark-probe {
  display: none;

  &[data-n="1"] {
    color: var(--msui-logo-1);
  }

  &[data-n="2"] {
    color: var(--msui-logo-2);
  }

  &[data-n="3"] {
    color: var(--msui-logo-3);
  }

  &[data-n="4"] {
    color: var(--msui-logo-4);
  }

  &[data-n="0"] {
    color: var(--msui-text-subtle);
  }
}
}
`);
