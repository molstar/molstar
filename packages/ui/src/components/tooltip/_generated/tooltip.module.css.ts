/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from tooltip.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'muted': 'msui-tooltip-muted',
  'tooltip': 'msui-tooltip-tooltip',
  'tooltip-in': 'msui-tooltip-tooltip-in',
} as const;

defineStyles('components/tooltip/tooltip.module.css', `@layer msui.components {
.msui-tooltip-tooltip {
  z-index: var(--msui-z-tooltip);
  max-width: 240px;
  padding: var(--msui-space-1) var(--msui-space-2);
  border-radius: var(--msui-radius-md);
  background: var(--msui-text);
  color: var(--msui-bg);
  font-size: var(--msui-text-xs);
  white-space: nowrap;
  text-overflow: ellipsis;
  pointer-events: none;
  animation: msui-tooltip-tooltip-in var(--msui-motion-fast) var(--msui-ease);
  line-height: 1.3;
  position: fixed;
  top: 0;
  left: 0;
  overflow: hidden;

  &[hidden] {
    display: none;
  }
}

.msui-tooltip-tooltip[data-kind="row"] {
  max-width: 480px;
  padding: var(--msui-space-1) var(--msui-control-inset);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  color: var(--msui-text);
  box-shadow: 0 0 0 1px var(--msui-line),
    var(--msui-shadow-float);
  font-size: var(--msui-text-sm);
  white-space: normal;
  overflow-wrap: anywhere;
  line-height: 1.4;
}

.msui-tooltip-muted {
  color: var(--msui-text-muted);
}

@keyframes msui-tooltip-tooltip-in {
  from {
    opacity: 0;
  }
}
}
`);
