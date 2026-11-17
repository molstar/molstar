/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from placeholders.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'band': 'msui-stage-placeholders-band',
  'boxes': 'msui-stage-placeholders-boxes',
  'canvas': 'msui-stage-placeholders-canvas',
  'column': 'msui-stage-placeholders-column',
  'letters': 'msui-stage-placeholders-letters',
  'residues': 'msui-stage-placeholders-residues',
  'ruler': 'msui-stage-placeholders-ruler',
  'selection': 'msui-stage-placeholders-selection',
  'structure': 'msui-stage-placeholders-structure',
} as const;

defineStyles('stage/placeholders.module.css', `@layer msui.components {
.msui-stage-placeholders-canvas {
  position: absolute;
  inset: 0;
}

.msui-stage-placeholders-structure {
  width: 100%;
  height: 100%;
  display: block;
}

.msui-stage-placeholders-band {
  padding: var(--msui-space-1) var(--msui-space-2) 0;
  scrollbar-width: thin;
  overflow: auto hidden;
}

.msui-stage-placeholders-residues {
  font-family: var(--msui-font-mono);
  font-size: var(--msui-text-xs);
  display: block;
}

.msui-stage-placeholders-ruler {
  fill: var(--msui-text-subtle);
  stroke: var(--msui-line-strong);
  text-anchor: middle;

  & text {
    stroke: none;
  }
}

.msui-stage-placeholders-boxes rect {
  fill: color-mix(in oklab, var(--msui-text) 8%, transparent);

  &[data-kind="hydrophobic"] {
    fill: color-mix(in oklab, var(--msui-warning) 22%, transparent);
  }

  &[data-kind="positive"] {
    fill: color-mix(in oklab, var(--msui-info) 26%, transparent);
  }

  &[data-kind="negative"] {
    fill: color-mix(in oklab, var(--msui-danger) 24%, transparent);
  }

  &[data-kind="polar"] {
    fill: color-mix(in oklab, var(--msui-success) 22%, transparent);
  }
}

.msui-stage-placeholders-letters {
  fill: var(--msui-text);
  dominant-baseline: central;
  text-anchor: middle;
  pointer-events: none;
}

.msui-stage-placeholders-selection {
  fill: none;
  stroke: var(--msui-accent);
  stroke-width: 2px;
}

.msui-stage-placeholders-column {
  fill: color-mix(in oklab, var(--msui-text) 12%, transparent);
  pointer-events: none;
  visibility: hidden;

  &[data-shown] {
    visibility: visible;
  }
}
}
`);
