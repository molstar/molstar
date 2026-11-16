/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from segmented.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'head': 'msui-segmented-head',
  'label': 'msui-segmented-label',
  'segment': 'msui-segmented-segment',
  'segmented': 'msui-segmented-segmented',
  'value': 'msui-segmented-value',
  'well': 'msui-segmented-well',
} as const;

defineStyles('components/segmented/segmented.module.css', `@layer msui.components {
.msui-segmented-segmented {
  --msui-segmented-height: var(--msui-row);
  --msui-segmented-pad: var(--msui-surface-pad);
  --msui-segmented-radius: var(--msui-radius-lg);
  --msui-segmented-icon: var(--msui-icon);
  gap: var(--msui-space-1);
  min-width: 0;
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  grid-template-columns: minmax(0, 1fr);
  display: grid;

  &[data-size="sm"] {
    --msui-segmented-height: var(--msui-control-height);
    --msui-segmented-pad: 2px;
    --msui-segmented-radius: var(--msui-radius-md);
    --msui-segmented-icon: var(--msui-icon-sm);
  }
}

.msui-segmented-head {
  column-gap: var(--msui-space-2);
  flex-wrap: wrap;
  min-width: 0;
  display: flex;
}

.msui-segmented-label {
  color: var(--msui-text-muted);
  white-space: nowrap;
}

.msui-segmented-value {
  color: var(--msui-text);
  white-space: nowrap;
}

.msui-segmented-well {
  padding: var(--msui-segmented-pad);
  border-radius: var(--msui-segmented-radius);
  background: var(--msui-surface-sunken);
  box-shadow: var(--msui-surface-edge);
  grid-auto-columns: minmax(max-content, 1fr);
  grid-auto-flow: column;
  display: grid;
}

.msui-segmented-segment {
  --msui-icon-size: var(--msui-segmented-icon);
  justify-content: center;
  align-items: center;
  gap: var(--msui-icon-gap);
  min-width: var(--msui-segmented-height);
  height: var(--msui-segmented-height);
  padding: 0 var(--msui-control-inset);
  border-radius: max(0px, calc(var(--msui-segmented-radius) - var(--msui-segmented-pad)));
  color: var(--msui-text-muted);
  white-space: nowrap;
  display: flex;

  [data-icons] > & {
    padding: 0;
  }

  &:hover:not([aria-disabled="true"]) {
    color: var(--msui-text);
  }

  &[aria-checked="true"] {
    background: var(--msui-surface-raised);
    box-shadow: 0 0 0 1px var(--msui-line),
      inset 0 0 0 1px var(--msui-selected-edge);
    color: var(--msui-accent-strong);
  }

  &[aria-disabled="true"] {
    color: color-mix(in oklab, var(--msui-text-subtle), transparent 50%);
    cursor: default;
  }
}
}
`);
