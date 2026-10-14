/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sequence-control.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'root': 'msp-viewer-sequence-control-root',
  'surface': 'msp-viewer-sequence-control-surface',
} as const;

defineStyles('viewer/sequence/sequence-control.module.css', `@layer msp.components {
.msp-viewer-sequence-control-root {
  width: 100%;
  height: var(--msp-dock-size, 100%);
  box-sizing: border-box;
  min-height: 0;
  padding-top: var(--msp-space-3);
  padding-bottom: var(--msp-band-fade, 0px);
  color: var(--msp-color-text);
  pointer-events: none;
  flex-direction: column;
  display: flex;
  position: relative;
}

.msp-viewer-sequence-control-surface {
  contain: layout paint;
  pointer-events: auto;
  flex: 1 1 0;
  min-height: 0;

  &[data-fade] {
    mask-image: linear-gradient(to bottom,
      transparent,
      #000 var(--seq-fade-top, 0px),
      #000 calc(100% - var(--seq-fade-bottom, 0px)),
      transparent);
  }
}

.msp-viewer-sequence-control-surface [data-pass^="vnav"] {
  opacity: 0;
  transition: opacity var(--msp-motion-base) var(--msp-ease);
}

.msp-viewer-sequence-control-surface[data-bars] [data-pass^="vnav"] {
  opacity: 1;
}
}
`);
