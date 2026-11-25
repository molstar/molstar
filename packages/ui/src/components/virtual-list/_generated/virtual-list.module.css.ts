/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from virtual-list.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'content': 'msui-virtual-list-content',
  'enter': 'msui-virtual-list-enter',
  'glide': 'msui-virtual-list-glide',
  'list': 'msui-virtual-list-list',
  'row': 'msui-virtual-list-row',
  'window': 'msui-virtual-list-window',
} as const;

defineStyles('components/virtual-list/virtual-list.module.css', `@layer msui.components {
.msui-virtual-list-list {
  flex: 1 1 0;
}

.msui-virtual-list-content {
  contain: strict;
  background-image: linear-gradient(to bottom,
    transparent calc(50% - 4px),
    color-mix(in oklab, var(--msui-text) 7%, transparent) calc(50% - 4px),
    color-mix(in oklab, var(--msui-text) 7%, transparent) calc(50% + 4px),
    transparent calc(50% + 4px));
  background-size: calc(100% - 2 * var(--msui-inner) - var(--msui-indent)) var(--msui-virtual-row, var(--msui-row));
  background-repeat: repeat-y;
  background-position: calc(var(--msui-inner) + var(--msui-indent)) 0;
  width: 100%;
  position: relative;
}

.msui-virtual-list-content[data-revealing] {
  background-image: none;
}

.msui-virtual-list-window {
  will-change: transform;
  background-color: var(--msui-virtual-bg, var(--msui-surface));
  position: absolute;
  inset: 0 0 auto;
}

.msui-virtual-list-row {
  contain: layout style;
  background-color: var(--msui-virtual-bg, var(--msui-surface));
  position: absolute;
  inset: 0 0 auto;
}

.msui-virtual-list-glide {
  z-index: 1;
  transition: transform var(--msui-motion-base) var(--msui-ease-in-out);
}

.msui-virtual-list-enter {
  animation: msui-reveal-fade var(--msui-motion-base) var(--msui-ease-in-out);
}

.msui-virtual-list-list[data-edge] .msui-virtual-list-window, .msui-virtual-list-list[data-edge] .msui-virtual-list-row {
  right: var(--msui-scroll-edge, var(--msui-space-1));
}
}
`);
