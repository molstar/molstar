/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from expandable.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'body': 'msui-expandable-body',
  'chevron': 'msui-expandable-chevron',
  'expandable': 'msui-expandable-expandable',
  'head': 'msui-expandable-head',
  'header': 'msui-expandable-header',
  'list': 'msui-expandable-list',
  'title': 'msui-expandable-title',
  'trailing': 'msui-expandable-trailing',
} as const;

defineStyles('components/expandable/expandable.module.css', `@layer msui.components {
.msui-expandable-list {
  gap: var(--msui-space-1);
  min-width: 0;
  display: grid;
}

.msui-expandable-expandable {
  border-radius: var(--msui-radius-md);
  min-width: 0;
  display: grid;

  &[data-variant="well"][data-open] {
    background: var(--msui-surface-sunken);
  }
}

.msui-expandable-head {
  display: contents;
}

.msui-expandable-header {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: center;
  gap: var(--msui-space-2);
  width: 100%;
  min-width: 0;
  height: var(--msui-header);
  padding-inline: var(--msui-inner) calc(var(--msui-inner) + var(--msui-space-1));
  border-radius: var(--msui-radius-md);
  color: var(--msui-text);
  font-size: var(--msui-text-md);
  font-weight: var(--msui-weight-strong);
  line-height: var(--msui-leading-md);
  text-align: start;
  transition: background-color var(--msui-motion-fast) var(--msui-ease);
  display: flex;

  &:hover {
    background: var(--msui-control-hover);
  }

  [data-open] > .msui-expandable-head > & {
    border-end-end-radius: 0;
    border-end-start-radius: 0;
  }
}

.msui-expandable-title {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  flex: 1;
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-expandable-trailing {
  color: var(--msui-text-muted);
  font-size: var(--msui-text-sm);
  font-weight: var(--msui-weight-regular);
  line-height: var(--msui-leading-sm);
  flex: none;
}

.msui-expandable-chevron {
  color: var(--msui-text-subtle);
  transition: transform var(--msui-motion-fast) var(--msui-ease);
  flex: none;

  [aria-expanded="false"] > & {
    transform: rotate(-90deg);
  }

  [dir="rtl"] [aria-expanded="false"] > & {
    transform: rotate(90deg);
  }
}

.msui-expandable-body {
  min-width: 0;
  padding-block: var(--msui-space-1) var(--msui-space-2);

  &[hidden] {
    display: none;
  }
}

.msui-expandable-expandable[data-variant="flat"] {
  & .msui-expandable-header {
    height: var(--msui-row);
    font-size: var(--msui-text-sm);
    line-height: var(--msui-leading-sm);

    &:hover {
      color: var(--msui-accent);
      background: none;
    }
  }

  & .msui-expandable-body {
    padding-block: 0 var(--msui-space-2);
  }
}
}
`);
