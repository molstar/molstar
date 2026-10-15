/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from track-gutter.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'gutter': 'msp-viewer-track-gutter-gutter',
  'icon': 'msp-viewer-track-gutter-icon',
} as const;

defineStyles('viewer/sequence/track-gutter.module.css', `@layer msp.components {
.msp-viewer-track-gutter-gutter {
  top: var(--msp-space-3);
  pointer-events: none;
  width: 28px;
  position: absolute;
  left: 0;
}

.msp-viewer-track-gutter-icon {
  width: 24px;
  min-height: max(24px, var(--msp-touch-target, 0px));
  margin-left: var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  opacity: 0;
  pointer-events: none;
  justify-content: center;
  align-items: center;
  display: flex;
  position: absolute;
  top: 0;
  left: 0;

  &[data-hover], &[data-active] {
    opacity: 1;
    pointer-events: auto;
  }

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[data-active] {
    color: var(--msp-color-accent);
  }
}
}
`);
