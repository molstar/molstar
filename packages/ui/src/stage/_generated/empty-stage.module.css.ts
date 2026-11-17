/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from empty-stage.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'empty': 'msui-stage-empty-stage-empty',
} as const;

defineStyles('stage/empty-stage.module.css', `@layer msui.components {
.msui-stage-empty-stage-empty {
  inset: var(--msui-stage-top-inset, 0px) 0 0;
  pointer-events: none;
  opacity: 1;
  transition: opacity var(--msui-motion-slow) var(--msui-ease);
  place-items: center;
  display: grid;
  position: absolute;

  @starting-style {
    opacity: 0;
  }

  & > * {
    pointer-events: auto;
  }

  &[data-phase="out"] {
    opacity: 0;
    transition: opacity var(--msui-stage-mark-out, .25s) var(--msui-ease-in);

    & > * {
      pointer-events: none;
    }
  }

  @media (prefers-reduced-motion: reduce) {
    transition: none;

    &[data-phase="out"] {
      transition: none;
    }
  }
}
}
`);
