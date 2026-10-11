/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from scene-tree.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'action': 'msp-structure-scene-tree-action',
  'actions': 'msp-structure-scene-tree-actions',
  'description': 'msp-structure-scene-tree-description',
  'empty': 'msp-structure-scene-tree-empty',
  'expander': 'msp-structure-scene-tree-expander',
  'hidden': 'msp-structure-scene-tree-hidden',
  'keep': 'msp-structure-scene-tree-keep',
  'kind': 'msp-structure-scene-tree-kind',
  'label': 'msp-structure-scene-tree-label',
  'leaf': 'msp-structure-scene-tree-leaf',
  'menu': 'msp-structure-scene-tree-menu',
  'menuItem': 'msp-structure-scene-tree-menuItem',
  'picked': 'msp-structure-scene-tree-picked',
  'row': 'msp-structure-scene-tree-row',
  'subChevron': 'msp-structure-scene-tree-subChevron',
  'subItems': 'msp-structure-scene-tree-subItems',
  'subMenu': 'msp-structure-scene-tree-subMenu',
  'tree': 'msp-structure-scene-tree-tree',
} as const;

defineStyles('views/structure/scene-tree.module.css', `@layer msp.components {
.msp-structure-scene-tree-tree {
  display: grid;
}

.msp-structure-scene-tree-empty {
  padding: var(--msp-space-4);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
}

.msp-structure-scene-tree-row {
  align-items: center;
  gap: var(--msp-space-2);
  height: 26px;
  padding: 0 var(--msp-space-1) 0 calc(var(--msp-space-1) + var(--msp-tree-depth, 0) * 14px);
  border-radius: var(--msp-radius-md);
  cursor: default;
  user-select: none;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-structure-scene-tree-picked {
  background: var(--msp-color-control);

  & .msp-structure-scene-tree-label {
    color: var(--msp-color-accent);
  }
}

.msp-structure-scene-tree-expander {
  border-radius: var(--msp-radius-sm);
  width: 18px;
  height: 18px;
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  place-items: center;
  display: grid;

  [aria-expanded="true"] > .msp-structure-scene-tree-row > & {
    transform: rotate(90deg);
  }

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-structure-scene-tree-leaf {
  visibility: hidden;
}

.msp-structure-scene-tree-kind {
  color: var(--msp-color-text-muted);
}

.msp-structure-scene-tree-label {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: none;
  max-width: 60%;
  font-weight: 500;
  overflow: hidden;
}

.msp-structure-scene-tree-description {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  overflow: hidden;
}

.msp-structure-scene-tree-actions {
  gap: 1px;
  display: flex;
}

.msp-structure-scene-tree-action {
  border-radius: var(--msp-radius-sm);
  width: 22px;
  height: 22px;
  color: var(--msp-color-text-muted);
  opacity: 0;
  place-items: center;
  display: grid;

  .msp-structure-scene-tree-row:hover &, &:focus-visible {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-structure-scene-tree-keep {
  opacity: 1;
}

.msp-structure-scene-tree-hidden {
  & .msp-structure-scene-tree-label, & .msp-structure-scene-tree-kind {
    color: var(--msp-color-text-subtle);
  }
}

.msp-structure-scene-tree-menu {
  margin: 2px var(--msp-space-1) var(--msp-space-1) calc(var(--msp-space-1) + var(--msp-tree-depth, 0) * 14px + 22px);
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  display: grid;
}

.msp-structure-scene-tree-menuItem {
  align-items: center;
  gap: var(--msp-space-3);
  height: 26px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-muted);
  text-align: left;
  grid-template-columns: 14px 1fr;
  display: grid;

  &:hover, &:focus-visible {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-structure-scene-tree-subChevron {
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  justify-self: center;

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}

.msp-structure-scene-tree-subItems {
  padding-left: calc(14px + var(--msp-space-3));
  display: grid;
}

.msp-structure-scene-tree-subMenu {
  display: grid;
}
}
`);
