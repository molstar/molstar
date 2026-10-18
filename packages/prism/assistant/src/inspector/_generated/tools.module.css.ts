/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from tools.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'field': 'msp-assistant-tools-field',
  'fields': 'msp-assistant-tools-fields',
  'toolList': 'msp-assistant-tools-toolList',
} as const;

defineStyles('assistant/inspector/tools.module.css', `@layer msp.components {
.msp-assistant-tools-toolList {
  gap: var(--msp-space-2);
  font-size: var(--msp-font-size-xs);
  grid-template-columns: minmax(0, 1fr);
  margin: 0;
  padding: 0;
  line-height: 1.45;
  list-style: none;
  display: grid;

  & li {
    align-items: baseline;
    gap: 0 var(--msp-space-3);
    flex-wrap: wrap;
    display: flex;
  }
}

.msp-assistant-tools-fields {
  gap: var(--msp-space-4);
  margin: 0;
  display: grid;
}

.msp-assistant-tools-field {
  gap: var(--msp-space-1);
  line-height: 1.45;
  display: grid;

  & dt {
    overflow-wrap: anywhere;
  }

  & dd {
    overflow-wrap: anywhere;
    margin: 0;
  }
}
}
`);
