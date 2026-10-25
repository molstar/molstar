/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from look-popover.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'button': 'msp-viewer-look-popover-button',
  'dot': 'msp-viewer-look-popover-dot',
  'mark': 'msp-viewer-look-popover-mark',
  'override': 'msp-viewer-look-popover-override',
  'overrideCount': 'msp-viewer-look-popover-overrideCount',
  'overrideLabel': 'msp-viewer-look-popover-overrideLabel',
  'overrides': 'msp-viewer-look-popover-overrides',
  'overridesTitle': 'msp-viewer-look-popover-overridesTitle',
} as const;

defineStyles('viewer/views/structure/look-popover.module.css', `@layer msp.components {
.msp-viewer-look-popover-button {
  position: relative;
}

.msp-viewer-look-popover-dot {
  background: var(--msp-color-accent);
  width: 6px;
  height: 6px;
  box-shadow: 0 0 0 1.5px var(--msp-color-raised);
  pointer-events: none;
  border-radius: 50%;
  position: absolute;
  top: 3px;
  right: 3px;
}

.msp-viewer-look-popover-overrides {
  gap: var(--msp-space-1);
  padding-top: var(--msp-space-2);
  border-top: 1px solid var(--msp-color-border);
  display: grid;
}

.msp-viewer-look-popover-overridesTitle {
  align-items: center;
  gap: var(--msp-space-1);
  color: var(--msp-color-text-subtle);
  display: flex;
}

.msp-viewer-look-popover-override {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 24px;
  display: flex;
}

.msp-viewer-look-popover-overrideLabel {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: auto;
  min-width: 0;
  overflow: hidden;
}

.msp-viewer-look-popover-overrideCount {
  color: var(--msp-color-text-subtle);
}

.msp-viewer-look-popover-mark {
  color: var(--msp-color-accent);
  flex: none;
}
}
`);
