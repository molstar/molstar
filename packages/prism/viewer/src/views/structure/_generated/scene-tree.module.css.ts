/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from scene-tree.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'action': 'msp-viewer-scene-tree-action',
  'actionSlot': 'msp-viewer-scene-tree-actionSlot',
  'actions': 'msp-viewer-scene-tree-actions',
  'description': 'msp-viewer-scene-tree-description',
  'empty': 'msp-viewer-scene-tree-empty',
  'expander': 'msp-viewer-scene-tree-expander',
  'hidden': 'msp-viewer-scene-tree-hidden',
  'inline': 'msp-viewer-scene-tree-inline',
  'inlineIcons': 'msp-viewer-scene-tree-inlineIcons',
  'inlineItem': 'msp-viewer-scene-tree-inlineItem',
  'item': 'msp-viewer-scene-tree-item',
  'keep': 'msp-viewer-scene-tree-keep',
  'kind': 'msp-viewer-scene-tree-kind',
  'label': 'msp-viewer-scene-tree-label',
  'leaf': 'msp-viewer-scene-tree-leaf',
  'loading': 'msp-viewer-scene-tree-loading',
  'look': 'msp-viewer-scene-tree-look',
  'lookInside': 'msp-viewer-scene-tree-lookInside',
  'menu': 'msp-viewer-scene-tree-menu',
  'menuItem': 'msp-viewer-scene-tree-menuItem',
  'more': 'msp-viewer-scene-tree-more',
  'moreCount': 'msp-viewer-scene-tree-moreCount',
  'picked': 'msp-viewer-scene-tree-picked',
  'row': 'msp-viewer-scene-tree-row',
  'rowBody': 'msp-viewer-scene-tree-rowBody',
  'short': 'msp-viewer-scene-tree-short',
  'subChevron': 'msp-viewer-scene-tree-subChevron',
  'subItems': 'msp-viewer-scene-tree-subItems',
  'subMenu': 'msp-viewer-scene-tree-subMenu',
  'swatch': 'msp-viewer-scene-tree-swatch',
  'tree': 'msp-viewer-scene-tree-tree',
} as const;

defineStyles('viewer/views/structure/scene-tree.module.css', `@layer msp.components {
.msp-viewer-scene-tree-tree {
  --msp-tree-indent: 14px;
}

.msp-viewer-scene-tree-empty {
  padding: var(--msp-space-4);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
}

.msp-viewer-scene-tree-item {
  outline: none;
  min-width: 0;
}

.msp-viewer-scene-tree-item:focus-visible > .msp-viewer-scene-tree-row, .msp-viewer-scene-tree-item:focus-visible > .msp-viewer-scene-tree-more, .msp-viewer-scene-tree-item:focus-visible > .msp-viewer-scene-tree-loading {
  box-shadow: inset 0 0 0 2px var(--msp-color-focus-ring);
}

.msp-viewer-scene-tree-row {
  align-items: center;
  gap: var(--msp-space-2);
  min-width: 0;
  height: max(26px, var(--msp-touch-target, 0px));
  padding: 0 var(--msp-space-1) 0
    calc(var(--msp-space-1) + var(--msp-tree-depth, 0) * var(--msp-tree-indent, var(--msp-icon-sm)));
  border-radius: var(--msp-radius-md);
  cursor: default;
  user-select: none;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-scene-tree-picked {
  background: var(--msp-color-control);

  & .msp-viewer-scene-tree-label {
    color: var(--msp-color-accent);
  }
}

.msp-viewer-scene-tree-expander {
  width: max(24px, var(--msp-touch-target, 0px));
  height: max(24px, var(--msp-touch-target, 0px));
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  flex: none;
  place-items: center;
  display: grid;

  [aria-expanded="true"] > .msp-viewer-scene-tree-row > & {
    transform: rotate(90deg);
  }

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-viewer-scene-tree-leaf {
  visibility: hidden;
}

.msp-viewer-scene-tree-kind {
  color: var(--msp-color-text-muted);
  flex: none;
}

.msp-viewer-scene-tree-label {
  min-width: 3ch;
  max-width: 60%;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 0 auto;
  font-weight: 500;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-scene-tree-description {
  min-width: 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-scene-tree-look {
  color: var(--msp-color-accent);
  flex: none;
  place-items: center;
  display: inline-grid;
}

.msp-viewer-scene-tree-lookInside {
  opacity: .5;
}

.msp-viewer-scene-tree-swatch {
  background: conic-gradient(var(--msp-accent-red) 0 25%,
    var(--msp-accent-orange) 0 50%,
    var(--msp-accent-green) 0 75%,
    var(--msp-accent-blue) 0);
  border-radius: 50%;
  flex: none;
  width: 8px;
  height: 8px;
}

.msp-viewer-scene-tree-actions {
  flex: none;
  gap: 1px;
  display: flex;
}

.msp-viewer-scene-tree-action {
  width: max(24px, var(--msp-touch-target, 0px));
  height: max(24px, var(--msp-touch-target, 0px));
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-muted);
  opacity: 0;
  place-items: center;
  display: grid;

  @media (hover: none) {
    opacity: 1;
  }

  .msp-viewer-scene-tree-row:hover &, &:focus-visible {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-viewer-scene-tree-actionSlot {
  width: max(24px, var(--msp-touch-target, 0px));
  place-items: center;
  display: grid;
}

.msp-viewer-scene-tree-keep {
  opacity: 1;
}

.msp-viewer-scene-tree-hidden {
  & .msp-viewer-scene-tree-label, & .msp-viewer-scene-tree-kind {
    color: var(--msp-color-text-subtle);
  }
}

.msp-viewer-scene-tree-menu {
  margin: var(--msp-space-1) var(--msp-space-1) var(--msp-space-1)
    calc(var(--msp-space-1) + var(--msp-tree-depth, 0) * var(--msp-tree-indent, var(--msp-icon-sm)) + 22px);
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  display: grid;

  &[data-popover], &[data-flat] {
    background: none;
    margin: 0;
    padding: 0;
  }
}

.msp-viewer-scene-tree-menuItem {
  align-items: center;
  gap: var(--msp-space-3);
  height: max(26px, var(--msp-touch-target, 0px));
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

.msp-viewer-scene-tree-subChevron {
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  justify-self: center;

  [aria-expanded="true"] > & {
    transform: rotate(90deg);
  }
}

.msp-viewer-scene-tree-subItems {
  padding-left: calc(var(--msp-icon-sm) + var(--msp-space-3));
  display: grid;
}

.msp-viewer-scene-tree-subMenu {
  display: grid;
}

.msp-viewer-scene-tree-more {
  align-items: baseline;
  gap: var(--msp-space-2);
  width: 100%;
  height: max(26px, var(--msp-touch-target, 0px));
  padding-left: calc(var(--msp-space-1) +
    var(--msp-tree-depth, 0) *
    var(--msp-tree-indent, var(--msp-icon-sm)) +
    22px);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-accent);
  font-size: var(--msp-font-size-sm);
  line-height: max(26px, var(--msp-touch-target, 0px));
  text-align: left;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-viewer-scene-tree-moreCount {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-viewer-scene-tree-loading {
  height: max(26px, var(--msp-touch-target, 0px));
  padding-left: calc(var(--msp-space-1) +
    var(--msp-tree-depth, 0) *
    var(--msp-tree-indent, var(--msp-icon-sm)) +
    40px);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  line-height: max(26px, var(--msp-touch-target, 0px));
}

.msp-viewer-scene-tree-inlineItem {
  cursor: default;
  grid-template-columns: 14px 1fr auto;

  &:hover {
    background: none;
  }
}

.msp-viewer-scene-tree-inlineIcons {
  gap: 2px;
  display: flex;
}

.msp-viewer-scene-tree-menuItem[data-danger] {
  color: var(--msp-color-danger);

  &:not(:first-child) {
    margin-top: var(--msp-space-1);
    box-shadow: 0 -1px 0 var(--msp-glass-edge);
    border-radius: 0 0 var(--msp-radius-sm) var(--msp-radius-sm);
  }

  &:hover, &:focus-visible {
    background: color-mix(in oklab, var(--msp-color-danger), transparent 88%);
    color: var(--msp-color-danger);
  }
}

.msp-viewer-scene-tree-rowBody {
  gap: var(--msp-space-2);
  min-width: 0;
  padding: var(--msp-space-1) var(--msp-space-2) var(--msp-space-2)
    calc(var(--msp-space-1) + var(--msp-tree-depth, 0) * var(--msp-tree-indent, var(--msp-icon-sm)) + 22px);
  display: grid;
}

.msp-viewer-scene-tree-description.msp-viewer-scene-tree-short {
  flex: none;
}

.msp-viewer-scene-tree-inline {
  min-width: 0;
  flex: 1;
  align-items: center;
  gap: 2px;
  margin-inline-start: var(--msp-space-2);
  display: flex;
}
}
`);
