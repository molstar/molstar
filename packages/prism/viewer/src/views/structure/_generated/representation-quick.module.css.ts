/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from representation-quick.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'quick': 'msp-viewer-representation-quick-quick',
} as const;

defineStyles('viewer/views/structure/representation-quick.module.css', `@layer msp.components {
.msp-viewer-representation-quick-quick {
  gap: var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}
}
`);
