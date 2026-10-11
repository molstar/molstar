/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from settings-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'body': 'msp-settings-settings-tab-body',
  'chevron': 'msp-settings-settings-tab-chevron',
  'group': 'msp-settings-settings-tab-group',
  'groups': 'msp-settings-settings-tab-groups',
  'header': 'msp-settings-settings-tab-header',
  'none': 'msp-settings-settings-tab-none',
  'search': 'msp-settings-settings-tab-search',
  'searchIcon': 'msp-settings-settings-tab-searchIcon',
  'searchInput': 'msp-settings-settings-tab-searchInput',
  'tab': 'msp-settings-settings-tab-tab',
  'title': 'msp-settings-settings-tab-title',
  'toggle': 'msp-settings-settings-tab-toggle',
} as const;

defineStyles('views/settings/settings-tab.module.css', `@layer msp.components {
.msp-settings-settings-tab-tab {
  height: 100%;
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  display: flex;
}

.msp-settings-settings-tab-search {
  align-items: center;
  gap: var(--msp-space-2);
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  flex: none;
  display: flex;
  position: relative;
}

.msp-settings-settings-tab-searchIcon {
  left: calc(var(--msp-space-2) + var(--msp-space-3));
  color: var(--msp-color-text-subtle);
  pointer-events: none;
  position: absolute;
}

.msp-settings-settings-tab-searchInput {
  min-width: 0;
  padding-left: calc(var(--msp-space-3) + 14px + var(--msp-space-2));
  flex: 1;
}

.msp-settings-settings-tab-groups {
  align-content: start;
  gap: var(--msp-space-1);
  min-height: 0;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  mask-image: linear-gradient(to bottom,
    transparent,
    #000 var(--msp-fade-top),
    #000 calc(100% - var(--msp-fade-bottom)),
    transparent);
  flex: 1 1 0;
  animation-name: msp-fade-scroll-top, msp-fade-scroll-bottom;
  animation-timing-function: linear;
  animation-fill-mode: both;
  animation-timeline: scroll(self);
  animation-range: 0 24px, calc(100% - 24px) 100%;
  display: grid;
  overflow: auto;
}

.msp-settings-settings-tab-group {
  border-radius: var(--msp-radius-md);

  &:has(.msp-settings-settings-tab-toggle[aria-expanded="true"]) {
    background: var(--msp-color-well);
  }
}

.msp-settings-settings-tab-header {
  border-radius: var(--msp-radius-md);
  align-items: center;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-settings-settings-tab-toggle {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  height: 32px;
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-muted);
  text-align: left;
  flex: 1;
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-settings-settings-tab-title {
  color: var(--msp-color-text);
  flex: 1;
  font-weight: 600;
}

.msp-settings-settings-tab-chevron {
  color: var(--msp-color-text-subtle);
  transition: transform var(--msp-motion-fast) var(--msp-ease);

  [aria-expanded="false"] > & {
    transform: rotate(-90deg);
  }
}

.msp-settings-settings-tab-body {
  --msp-label-width: 46%;
  padding: var(--msp-space-1) var(--msp-space-3) var(--msp-space-3);

  &[data-wide-labels] {
    --msp-label-width: calc(100% - 64px);
  }
}

.msp-settings-settings-tab-none {
  padding: var(--msp-space-4) var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  text-align: center;
  margin: 0;
}
}
`);
