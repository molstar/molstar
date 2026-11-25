/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from stage.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'band': 'msui-stage-band',
  'bottom': 'msui-stage-bottom',
  'canvas': 'msui-stage-canvas',
  'chips': 'msui-stage-chips',
  'corner': 'msui-stage-corner',
  'ink': 'msui-stage-ink',
  'lead': 'msui-stage-lead',
  'overlays': 'msui-stage-overlays',
  'slot': 'msui-stage-slot',
  'stage': 'msui-stage-stage',
  'status': 'msui-stage-status',
  'tools': 'msui-stage-tools',
  'top': 'msui-stage-top',
} as const;

defineStyles('stage/stage.module.css', `@layer msui.components {
.msui-stage-stage {
  --msui-stage-edge: var(--msui-space-2);
  --msui-stage-band-fade: var(--msui-space-5);
  --msui-stage-top-inset: 0px;
  background: var(--msui-stage-color, #fff);
  isolation: isolate;
  flex: auto;
  width: 100%;
  min-width: 0;
  height: 100%;
  min-height: 0;
  position: relative;
  overflow: hidden;
}

.msui-stage-canvas {
  opacity: 1;
  transition: opacity var(--msui-motion-slow) var(--msui-ease);
  position: absolute;
  inset: 0;

  &[data-quiet] {
    opacity: 0;
  }

  @media (prefers-reduced-motion: reduce) {
    transition: none;
  }
}

.msui-stage-ink {
  color: var(--msui-text);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  pointer-events: none;
  position: absolute;
  inset: 0;
}

.msui-stage-slot {
  pointer-events: none;
}

:where(.msui-stage-slot) > * {
  pointer-events: auto;
}

.msui-stage-top {
  pointer-events: none;
  position: absolute;
  top: 0;
  left: 0;
  right: 0;

  &:after {
    content: "";
    height: var(--msui-stage-band-fade);
    background: linear-gradient(to bottom, var(--msui-stage-color, #fff), transparent);
    display: block;
  }

  &[hidden] {
    display: none;
  }
}

.msui-stage-band {
  background: var(--msui-stage-color, #fff);
  pointer-events: auto;
}

.msui-stage-lead {
  top: calc(var(--msui-stage-top-inset) + var(--msui-stage-edge));
  left: var(--msui-stage-edge);
  justify-items: start;
  gap: var(--msui-space-2);
  max-width: calc(100% - 2 * var(--msui-stage-edge));
  display: grid;
  position: absolute;
}

.msui-stage-bottom {
  right: calc(var(--msui-stage-edge) + var(--msui-stage-inset-right, 0px));
  bottom: var(--msui-stage-edge);
  left: calc(var(--msui-stage-edge) + var(--msui-stage-inset-left, 0px));
  align-items: end;
  gap: var(--msui-space-2);
  pointer-events: none;
  grid-template-columns: minmax(0, 1fr) auto minmax(0, 1fr);
  grid-template-areas: "status chips corner";
  display: grid;
  position: absolute;
}

.msui-stage-status {
  align-items: end;
  gap: var(--msui-space-1);
  flex-wrap: wrap;
  grid-area: status;
  min-width: 0;
  display: flex;
}

.msui-stage-chips {
  justify-content: center;
  gap: var(--msui-space-1);
  grid-area: chips;
  display: flex;
}

.msui-stage-corner {
  --msui-hover-label-right: 0px;
  --msui-hover-label-bottom: calc(100% + var(--msui-space-2));
  grid-area: corner;
  justify-content: flex-end;
  align-items: end;
  min-width: 0;
  display: flex;
  position: relative;
}

.msui-stage-tools {
  transition: opacity var(--msui-motion-fast) var(--msui-ease);
  display: flex;

  .msui-stage-stage[data-toolbar="hidden"] & {
    opacity: 0;
    pointer-events: none;
  }

  .msui-stage-stage[data-toolbar="hidden"] &:focus-within {
    opacity: 1;
    pointer-events: auto;
  }
}

.msui-stage-stage[data-compact] {
  & .msui-stage-bottom {
    grid-template-columns: minmax(0, 1fr) auto;
    grid-template-areas: "status status"
                         "chips corner";
  }

  & .msui-stage-chips {
    justify-content: flex-start;
  }
}

.msui-stage-overlays {
  position: absolute;
  inset: 0;
}
}
`);
