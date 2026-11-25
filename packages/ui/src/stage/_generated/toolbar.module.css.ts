/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from toolbar.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'bar': 'msui-stage-toolbar-bar',
  'folded': 'msui-stage-toolbar-folded',
  'row': 'msui-stage-toolbar-row',
  'scroll': 'msui-stage-toolbar-scroll',
  'stageToolbar': 'msui-stage-toolbar-stageToolbar',
} as const;

defineStyles('stage/toolbar.module.css', `@layer msui.components {
.msui-stage-toolbar-stageToolbar {
  justify-content: flex-end;
  align-items: end;
  min-width: 0;
  max-width: 100%;
  display: flex;

  & [role="menu"] {
    max-height: var(--msui-stage-menu-room, calc(100vh - var(--msui-space-4)));
  }
}

.msui-stage-toolbar-bar {
  align-items: end;
  gap: var(--msui-space-2);
  min-width: 0;
  display: flex;

  &[hidden] {
    display: none;
  }
}

.msui-stage-toolbar-scroll {
  min-width: 0;
  margin: calc(-1 * var(--msui-space-3));
}

.msui-stage-toolbar-row {
  align-items: end;
  gap: var(--msui-space-2);
  width: max-content;
  padding: var(--msui-space-3);
  display: flex;
}

.msui-stage-toolbar-folded {
  display: flex;
}
}
`);
