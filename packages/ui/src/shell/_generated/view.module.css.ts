/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from view.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'actions': 'msui-shell-view-actions',
  'body': 'msui-shell-view-body',
  'bodyContent': 'msui-shell-view-bodyContent',
  'bodyHost': 'msui-shell-view-bodyHost',
  'bodyWrap': 'msui-shell-view-bodyWrap',
  'empty': 'msui-shell-view-empty',
  'emptyAction': 'msui-shell-view-emptyAction',
  'emptyText': 'msui-shell-view-emptyText',
  'error': 'msui-shell-view-error',
  'fill': 'msui-shell-view-fill',
  'fillHost': 'msui-shell-view-fillHost',
  'footer': 'msui-shell-view-footer',
  'footerSlot': 'msui-shell-view-footerSlot',
  'frame': 'msui-shell-view-frame',
  'header': 'msui-shell-view-header',
  'headerButton': 'msui-shell-view-headerButton',
  'heading': 'msui-shell-view-heading',
  'host': 'msui-shell-view-host',
  'stack': 'msui-shell-view-stack',
  'title': 'msui-shell-view-title',
  'toolbar': 'msui-shell-view-toolbar',
} as const;

defineStyles('shell/view.module.css', `@layer msui.components {
.msui-shell-view-frame {
  --msui-view-pad: 6px;
  width: 100%;
  min-width: 0;
  height: 100%;
  min-height: 0;
  color: var(--msui-text);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
  flex-direction: column;
  flex: auto;
  display: flex;
}

.msui-shell-view-header {
  align-items: center;
  gap: var(--msui-space-1);
  height: var(--msui-header);
  padding-inline: calc(var(--msui-view-pad) + var(--msui-inner));
  flex: none;
  min-width: 0;
  display: flex;

  &[data-head="rail"] {
    height: var(--msui-shell-head, var(--msui-header));
  }
}

.msui-shell-view-title {
  white-space: nowrap;
  min-width: 0;
  font-size: var(--msui-text-md);
  font-weight: var(--msui-weight-strong);
  line-height: var(--msui-leading-md);
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  flex: 1;
  margin: 0;
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-shell-view-actions {
  align-items: center;
  gap: var(--msui-space-1);
  flex: none;
  display: flex;
}

.msui-shell-view-headerButton {
  flex: none;
  display: inline-flex;
}

.msui-shell-view-bodyWrap {
  flex-direction: column;
  flex: auto;
  min-width: 0;
  min-height: 0;
  display: flex;
}

.msui-shell-view-body {
  flex: auto;
  min-height: 0;
}

.msui-shell-view-footerSlot {
  padding: var(--msui-space-1) var(--msui-view-pad);
  box-shadow: 0 -1px 0 var(--msui-shell-line, var(--msui-line));
  flex: none;

  &:empty {
    display: none;
  }
}

.msui-shell-view-bodyContent {
  min-height: 100%;
  padding: var(--msui-space-1) var(--msui-view-pad) var(--msui-view-pad);
  flex-direction: column;
  display: flex;
}

.msui-shell-view-frame:has( > .msui-shell-view-header [role="tablist"]) .msui-shell-view-bodyContent {
  padding-top: var(--msui-space-2);
}

.msui-shell-view-bodyHost {
  gap: var(--msui-space-2);
  flex-direction: column;
  flex: 1 0 auto;
  min-width: 0;
  display: flex;
}

.msui-shell-view-fill, .msui-shell-view-fillHost {
  flex-direction: column;
  flex: auto;
  min-height: 0;
  display: flex;
}

.msui-shell-view-host {
  min-width: 0;
}

.msui-shell-view-stack {
  gap: var(--msui-space-2);
  flex-direction: column;
  flex: 1 0 auto;
  min-width: 0;
  display: flex;

  &[data-gap="sections"] {
    gap: var(--msui-space-1);
  }
}

.msui-shell-view-toolbar {
  align-items: center;
  gap: var(--msui-space-1);
  min-height: var(--msui-row);
  padding-inline: var(--msui-inner);
  display: flex;
}

.msui-shell-view-empty {
  gap: var(--msui-space-2);
  padding: var(--msui-space-2) var(--msui-inner);
  justify-items: start;
  display: grid;
}

.msui-shell-view-emptyText {
  color: var(--msui-text-muted);
  margin: 0;
}

.msui-shell-view-emptyAction {
  gap: var(--msui-space-1);
  display: flex;
}

.msui-shell-view-error {
  --msui-icon-size: var(--msui-icon-sm);
  align-items: center;
  gap: var(--msui-space-1);
  min-height: var(--msui-row);
  padding-inline: var(--msui-inner);
  color: var(--msui-danger);
  display: flex;
}

.msui-shell-view-footer {
  min-width: 0;
}

.msui-shell-view-heading {
  padding: var(--msui-space-2) var(--msui-inner) 0;
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  font-weight: var(--msui-weight-strong);
  line-height: var(--msui-leading-xs);
  margin: 0;
}
}
`);
