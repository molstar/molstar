/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from dynamic-pdb.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'chip': 'msp-viewer-dynamic-pdb-chip',
} as const;

defineStyles('viewer/_stories/dynamic-pdb/dynamic-pdb.module.css', `@layer msp.components {
.msp-viewer-dynamic-pdb-chip {
  font-family: var(--msp-font-mono);
  font-size: var(--msp-font-size-xs);
  color: var(--msp-color-text);
  padding: 1px var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-well);
  cursor: pointer;
  border: none;
}

.msp-viewer-dynamic-pdb-chip:hover {
  background: var(--msp-color-hover);
}

.msp-viewer-dynamic-pdb-chip:focus-visible {
  outline: 2px solid var(--msp-color-focus-ring);
  outline-offset: 1px;
}

.msp-viewer-dynamic-pdb-chip:disabled {
  color: var(--msp-color-text-subtle);
  cursor: default;
}
}
`);
