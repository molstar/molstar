/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from skeleton.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'icon': 'msui-skeleton-icon',
  'in': 'msui-skeleton-in',
  'label': 'msui-skeleton-label',
  'layer': 'msui-skeleton-layer',
  'line': 'msui-skeleton-line',
  'loading': 'msui-skeleton-loading',
  'pulse': 'msui-skeleton-pulse',
  'round': 'msui-skeleton-round',
  'shape': 'msui-skeleton-shape',
  'swap': 'msui-skeleton-swap',
  'text': 'msui-skeleton-text',
} as const;

defineStyles('components/skeleton/skeleton.module.css', `@layer msui.components {
.msui-skeleton-loading {
  min-width: 0;
  animation: msui-skeleton-in var(--msui-motion-base) var(--msui-ease) var(--msui-skeleton-delay, .3s) backwards;
  display: grid;
}

@keyframes msui-skeleton-in {
  from {
    opacity: 0;
  }
}

.msui-skeleton-label {
  clip-path: inset(50%);
  white-space: nowrap;
  width: 1px;
  height: 1px;
  position: absolute;
  overflow: hidden;
}

.msui-skeleton-swap {
  min-width: 0;
  display: grid;

  & > * {
    grid-area: 1 / 1;
    min-width: 0;
  }
}

.msui-skeleton-layer {
  min-width: 0;
  display: grid;
}

.msui-skeleton-shape {
  border-radius: var(--msui-radius-sm);
  background: color-mix(in oklab, currentColor, transparent 90%);
  min-width: 0;
  animation: msui-skeleton-pulse 1.6s var(--msui-ease-in-out) var(--msui-skeleton-delay, .3s) infinite;
  display: block;
}

@keyframes msui-skeleton-pulse {
  50% {
    opacity: .5;
  }
}

@media (prefers-reduced-motion: reduce) {
  .msui-skeleton-shape {
    animation: none;
  }
}

.msui-skeleton-line {
  height: var(--msui-leading-sm);
  align-items: center;
  display: flex;

  & > .msui-skeleton-shape {
    height: var(--msui-space-2);
  }

  &[data-size="xs"] {
    height: var(--msui-leading-xs);
  }

  &[data-size="md"] {
    height: var(--msui-leading-md);

    & > .msui-skeleton-shape {
      height: var(--msui-space-3);
    }
  }
}

.msui-skeleton-icon {
  width: var(--msui-icon);
  height: var(--msui-icon);
  border-radius: var(--msui-radius-sm);
}

.msui-skeleton-text {
  height: var(--msui-space-2);
}

.msui-skeleton-round {
  border-radius: 50%;
}
}
`);
