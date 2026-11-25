/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from scroll-area.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'content': 'msui-scroll-area-content',
  'root': 'msui-scroll-area-root',
  'scrollbar': 'msui-scroll-area-scrollbar',
  'thumb': 'msui-scroll-area-thumb',
  'viewport': 'msui-scroll-area-viewport',
} as const;

defineStyles('components/scroll-area/scroll-area.module.css', `@layer msui.components {
.msui-scroll-area-root {
  flex-direction: column;
  min-width: 0;
  min-height: 0;
  display: flex;
  position: relative;
  overflow: hidden;
}

.msui-scroll-area-root[data-edge] {
  margin-inline-end: calc(-1 * var(--msui-scroll-edge, var(--msui-space-1)));
}

.msui-scroll-area-viewport {
  overscroll-behavior: contain;
  scrollbar-width: none;
  flex: auto;
  min-width: 0;
  min-height: 0;
  overflow: auto;

  &::-webkit-scrollbar {
    display: none;
  }

  .msui-scroll-area-root[data-chain] > & {
    overscroll-behavior: auto;
  }

  .msui-scroll-area-root[data-axis="y"] > & {
    overflow-x: hidden;
  }

  .msui-scroll-area-root[data-axis="x"] > & {
    overflow-y: hidden;
  }

  .msui-scroll-area-root:not([data-axis="x"], [data-no-fade]) > & {
    mask-image: linear-gradient(to bottom,
      transparent,
      #000 var(--msui-fade-top),
      #000 calc(100% - var(--msui-fade-bottom)),
      transparent);
    animation-name: msui-fade-scroll-top, msui-fade-scroll-bottom;
    animation-timing-function: linear;
    animation-fill-mode: both;
    animation-timeline: scroll(self);
    animation-range: 0 24px, calc(100% - 24px) 100%;
  }
}

.msui-scroll-area-content {
  min-width: 0;

  .msui-scroll-area-root[data-edge] > .msui-scroll-area-viewport > & {
    padding-inline-end: var(--msui-scroll-edge, var(--msui-space-1));
  }

  .msui-scroll-area-root[data-axis="x"] > .msui-scroll-area-viewport > &, .msui-scroll-area-root[data-axis="both"] > .msui-scroll-area-viewport > & {
    width: max-content;
    min-width: 100%;
  }
}

.msui-scroll-area-scrollbar {
  z-index: 1;
  opacity: 0;
  transition: opacity var(--msui-motion-base) var(--msui-ease);
  touch-action: none;
  user-select: none;
  display: flex;
  position: absolute;

  &[data-orientation="vertical"] {
    justify-content: center;
    width: 10px;
    top: 2px;
    bottom: 2px;
    right: 1px;
  }

  &[data-orientation="horizontal"] {
    flex-direction: column;
    justify-content: center;
    height: 10px;
    bottom: 1px;
    left: 2px;
    right: 2px;
  }

  .msui-scroll-area-root:hover > &, .msui-scroll-area-root:focus-within > &, .msui-scroll-area-root[data-scrolling] > &, &[data-dragging] {
    opacity: 1;
  }

  .msui-scroll-area-root[data-edge] > &[data-orientation="vertical"] {
    width: 8px;
    right: 0;
  }

  .msui-scroll-area-root:not([data-overflow-y]) > &[data-orientation="vertical"], .msui-scroll-area-root:not([data-overflow-x]) > &[data-orientation="horizontal"] {
    display: none;
  }
}

.msui-scroll-area-thumb {
  background: color-mix(in oklab, var(--msui-text), transparent 72%);
  transition: scale var(--msui-motion-fast) var(--msui-ease);
  border-radius: 3px;
  flex: none;

  [data-orientation="vertical"] > & {
    transform-origin: 100%;
    width: 4px;
  }

  [data-orientation="horizontal"] > & {
    transform-origin: bottom;
    height: 4px;
  }

  .msui-scroll-area-scrollbar:hover > &, [data-dragging] > & {
    background: color-mix(in oklab, var(--msui-text), transparent 52%);
  }

  [data-orientation="vertical"]:hover > &, [data-orientation="vertical"][data-dragging] > & {
    scale: 1.5 1;
  }

  [data-orientation="horizontal"]:hover > &, [data-orientation="horizontal"][data-dragging] > & {
    scale: 1 1.5;
  }
}
}
`);
