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
  'bar': 'msp-viewer-sequence-control-bar',
  'popup': 'msp-viewer-sequence-control-popup',
  'root': 'msp-viewer-sequence-control-root',
  'spacer': 'msp-viewer-sequence-control-spacer',
  'surface': 'msp-viewer-sequence-control-surface',
  'switch': 'msp-viewer-sequence-control-switch',
} as const;

defineStyles('viewer/sequence/sequence-control.module.css', `@layer msp.components {
.msp-viewer-sequence-control-root {
  width: 100%;
  height: var(--msp-dock-size, 100%);
  min-height: 0;
  color: var(--msp-color-text);
  flex-direction: column;
  display: flex;
  position: relative;
}

.msp-viewer-sequence-control-bar {
  align-items: center;
  gap: var(--msp-space-2);
  height: var(--msp-control-height);
  padding: 0 var(--msp-space-1);
  pointer-events: none;
  contain: layout paint;
  flex: none;
  display: flex;

  & > * {
    pointer-events: auto;
  }
}

.msp-viewer-sequence-control-switch {
  width: 140px;
}

.msp-viewer-sequence-control-spacer {
  pointer-events: none;
  flex: 1;
}

.msp-viewer-sequence-control-surface {
  contain: layout paint;
  flex: 1 1 0;
  min-height: 0;
}

.msp-viewer-sequence-control-popup {
  top: var(--msp-control-height);
  right: var(--msp-space-1);
  z-index: var(--msp-z-overlay);
  max-width: 280px;
  max-height: 320px;
  padding: var(--msp-space-3);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-panel-bg);
  box-shadow: var(--msp-glass-shadow);
  font-size: var(--msp-font-size-sm);
  position: absolute;
  overflow: auto;
}
}
`);
