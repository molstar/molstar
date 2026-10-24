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
  'said': 'msp-viewer-sequence-control-said',
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
  position: relative;

  &:not([data-pointer]):has( > :focus-visible):after {
    content: "";
    border-radius: var(--msp-radius-md);
    box-shadow: inset 0 0 0 2px var(--msp-color-focus-ring);
    pointer-events: none;
    position: absolute;
    inset: 0;
  }

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

.msp-viewer-sequence-control-said {
  clip-path: inset(50%);
  white-space: nowrap;
  width: 1px;
  height: 1px;
  position: absolute;
  overflow: hidden;
}
}
`);
