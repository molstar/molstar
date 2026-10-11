/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from remote-states-panel.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'action': 'msp-state-remote-states-panel-action',
  'caption': 'msp-state-remote-states-panel-caption',
  'captionRow': 'msp-state-remote-states-panel-captionRow',
  'date': 'msp-state-remote-states-panel-date',
  'empty': 'msp-state-remote-states-panel-empty',
  'item': 'msp-state-remote-states-panel-item',
  'list': 'msp-state-remote-states-panel-list',
  'menu': 'msp-state-remote-states-panel-menu',
  'menuItem': 'msp-state-remote-states-panel-menuItem',
  'name': 'msp-state-remote-states-panel-name',
  'open': 'msp-state-remote-states-panel-open',
  'panel': 'msp-state-remote-states-panel-panel',
  'row': 'msp-state-remote-states-panel-row',
  'server': 'msp-state-remote-states-panel-server',
  'spin': 'msp-state-remote-states-panel-spin',
  'states': 'msp-state-remote-states-panel-states',
  'upload': 'msp-state-remote-states-panel-upload',
} as const;

defineStyles('views/state/remote-states-panel.module.css', `@layer msp.components {
.msp-state-remote-states-panel-panel {
  gap: var(--msp-space-4);
  height: 100%;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  display: flex;
}

.msp-state-remote-states-panel-upload {
  flex: none;
}

.msp-state-remote-states-panel-caption {
  padding: var(--msp-space-1) 0;
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
  font-weight: 500;
}

.msp-state-remote-states-panel-states {
  gap: var(--msp-space-1);
  flex-direction: column;
  flex: 1 1 0;
  min-height: 0;
  display: flex;
}

.msp-state-remote-states-panel-captionRow {
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  justify-content: space-between;
  align-items: center;
  display: flex;
}

.msp-state-remote-states-panel-list {
  min-height: 0;
  mask-image: linear-gradient(to bottom,
    transparent,
    #000 var(--msp-fade-top),
    #000 calc(100% - var(--msp-fade-bottom)),
    transparent);
  flex: 1 1 0;
  grid-template-columns: minmax(0, 1fr);
  align-content: start;
  animation-name: msp-fade-scroll-top, msp-fade-scroll-bottom;
  animation-timing-function: linear;
  animation-fill-mode: both;
  animation-timeline: scroll(self);
  animation-range: 0 24px, calc(100% - 24px) 100%;
  display: grid;
  overflow: auto;
}

.msp-state-remote-states-panel-empty {
  padding: var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  margin: 0;
}

.msp-state-remote-states-panel-item {
  min-width: 0;
  display: grid;
}

.msp-state-remote-states-panel-row {
  border-radius: var(--msp-radius-md);
  align-items: center;
  min-width: 0;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }
}

.msp-state-remote-states-panel-open {
  align-items: center;
  gap: var(--msp-space-3);
  min-width: 0;
  height: 28px;
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-muted);
  text-align: left;
  flex: 1;
  display: flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-state-remote-states-panel-name {
  min-width: 3ch;
  color: var(--msp-color-text);
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  flex: 0 auto;
  font-weight: 500;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-state-remote-states-panel-date {
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

.msp-state-remote-states-panel-action {
  width: 22px;
  height: 22px;
  margin-right: var(--msp-space-1);
  border-radius: var(--msp-radius-sm);
  color: var(--msp-color-text-muted);
  opacity: 0;
  flex: none;
  place-items: center;
  display: grid;

  .msp-state-remote-states-panel-row:hover &, &:focus-visible, &[aria-expanded="true"] {
    opacity: 1;
  }

  &:hover {
    background: var(--msp-color-pressed);
    color: var(--msp-color-text);
  }
}

.msp-state-remote-states-panel-menu {
  margin: 2px var(--msp-space-1) var(--msp-space-1) calc(var(--msp-space-2) + 15px + var(--msp-space-3));
  padding: var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 40%);
  display: grid;
}

.msp-state-remote-states-panel-menuItem {
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

  &[data-danger] {
    color: var(--msp-color-danger);

    &:hover, &:focus-visible {
      background: color-mix(in oklab, var(--msp-color-danger), transparent 88%);
      color: var(--msp-color-danger);
    }
  }
}

.msp-state-remote-states-panel-server {
  align-items: center;
  gap: var(--msp-space-3);
  padding: 0 var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  flex: none;
  grid-template-columns: auto minmax(0, 1fr);
  display: grid;
}

.msp-state-remote-states-panel-spin {
  animation: .9s linear infinite msp-state-remote-states-panel-spin;
}

@keyframes msp-state-remote-states-panel-spin {
  to {
    transform: rotate(360deg);
  }
}
}
`);
