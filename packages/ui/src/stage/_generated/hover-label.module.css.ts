/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from hover-label.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'body': 'msui-stage-hover-label-body',
  'code': 'msui-stage-hover-label-code',
  'corner': 'msui-stage-hover-label-corner',
  'cursor': 'msui-stage-hover-label-cursor',
  'item': 'msui-stage-hover-label-item',
  'label': 'msui-stage-hover-label-label',
  'link': 'msui-stage-hover-label-link',
  'out': 'msui-stage-hover-label-out',
  'primary': 'msui-stage-hover-label-primary',
  'secondary': 'msui-stage-hover-label-secondary',
} as const;

defineStyles('stage/hover-label.module.css', `@layer msui.components {
.msui-stage-hover-label-label {
  z-index: var(--msui-z-tooltip);
  max-width: calc(70 * var(--msui-unit));
  padding: var(--msui-space-2) var(--msui-space-3);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-float);
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  pointer-events: none;
  transition: none;
  position: absolute;
}

.msui-stage-hover-label-cursor {
  transform: translate(var(--msui-hover-x, 0px), var(--msui-hover-y, 0px));
  will-change: transform;
  top: 0;
  left: 0;
}

.msui-stage-hover-label-corner {
  right: var(--msui-hover-label-right, var(--msui-space-3));
  bottom: var(--msui-hover-label-bottom, var(--msui-space-3));
  width: max-content;
}

.msui-stage-hover-label-body {
  --cap: var(--msui-hover-label-height, calc(64 * var(--msui-unit)));
  max-height: var(--cap);
  overflow-wrap: anywhere;
  mask-image: linear-gradient(to bottom, #000 calc(var(--cap) - 2.5em), transparent var(--cap));
  display: grid;
  overflow: hidden;
}

.msui-stage-hover-label-out {
  opacity: 0;
  visibility: hidden;
  transition: opacity var(--msui-motion-fast) var(--msui-ease-in),
    visibility 0s var(--msui-motion-fast);

  @media (prefers-reduced-motion: reduce) {
    transition: none;
  }
}

.msui-stage-hover-label-primary {
  font-weight: var(--msui-weight-strong);
}

.msui-stage-hover-label-secondary {
  color: var(--msui-text-muted);
}

.msui-stage-hover-label-item {
  text-indent: -1em;
  padding-inline-start: 1em;
}

.msui-stage-hover-label-link {
  color: var(--msui-accent-strong);
  text-decoration: underline;
}

.msui-stage-hover-label-code {
  font-family: var(--msui-font-mono);
  font-size: .95em;
}
}
`);
