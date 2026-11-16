/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from iso-level.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'inside': 'msui-iso-level-inside',
  'marker': 'msui-iso-level-marker',
  'off': 'msui-iso-level-off',
  'plot': 'msui-iso-level-plot',
  'track': 'msui-iso-level-track',
} as const;

defineStyles('components/iso-level/iso-level.module.css', `@layer msui.components {
.msui-iso-level-plot {
  cursor: ew-resize;
}

.msui-iso-level-inside {
  background: color-mix(in oklab, var(--msui-level, var(--msui-text-muted)), transparent 82%);
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msui-iso-level-track {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msui-iso-level-marker {
  top: 0;
  bottom: 0;
  left: calc(-1 * var(--msui-space-3) / 2);
  width: var(--msui-space-3);
  pointer-events: auto;
  cursor: ew-resize;
  outline: none;
  position: absolute;

  &:before {
    content: "";
    background: var(--msui-level, var(--msui-text));
    width: 2px;
    position: absolute;
    top: 0;
    bottom: 0;
    left: calc(50% - 1px);
  }

  &:after {
    content: "";
    top: var(--msui-space-1);
    left: calc(50% - var(--msui-space-1));
    width: var(--msui-space-2);
    height: var(--msui-space-2);
    background: var(--msui-level, var(--msui-text));
    box-shadow: 0 0 0 1px var(--msui-surface-sunken);
    transition: transform var(--msui-motion-fast) var(--msui-ease);
    border-radius: 50%;
    position: absolute;
  }

  &:hover:after, &:focus-visible:after {
    transform: scale(1.25);
  }

  &:focus-visible:after {
    box-shadow: 0 0 0 1px var(--msui-surface-sunken),
      0 0 0 2px var(--msui-focus);
  }
}

.msui-iso-level-off {
  opacity: .35;
}
}
`);
