/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from help.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'description': 'msui-params-help-description',
  'help': 'msui-params-help-help',
  'table': 'msui-params-help-table',
  'text': 'msui-params-help-text',
  'toggle': 'msui-params-help-toggle',
} as const;

defineStyles('components/params/help.module.css', `@layer msui.components {
.msui-params-help-text {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;

  &[data-help] {
    text-decoration: underline dotted color-mix(in oklab, currentColor, transparent 45%);
    text-underline-offset: 3px;
    cursor: help;
    text-decoration-thickness: 1px;
  }
}

.msui-params-help-toggle {
  --msui-icon-size: var(--msui-icon-sm);
  width: var(--msui-icon-slot);
  height: var(--msui-leading-sm);
  color: var(--msui-text-subtle);
  opacity: 0;
  border-radius: 50%;
  flex: none;
  justify-content: center;
  align-items: center;
  margin-inline-start: var(--msui-space-1);
  display: inline-flex;
  position: relative;

  &:after {
    content: "";
    width: max(var(--msui-row), var(--msui-touch-target, 0px));
    height: max(var(--msui-row), var(--msui-touch-target, 0px));
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
  }

  :where([data-field]:hover) &, &[aria-expanded="true"], &:focus-visible {
    opacity: 1;
  }

  &:hover, &[aria-expanded="true"] {
    color: var(--msui-text);
  }
}

.msui-params-help-help {
  grid-template-columns: subgrid;
  row-gap: var(--msui-space-2);
  margin-block: 0 var(--msui-space-1);
  margin-inline: calc(var(--msui-params-level, 0) * var(--msui-form-indent) - var(--msui-control-inset))
    calc(-1 * var(--msui-control-inset));
  padding: var(--msui-space-1) var(--msui-control-inset);
  border-radius: var(--msui-radius-md);
  background: color-mix(in oklab, var(--msui-control), transparent 40%);
  display: grid;
}

.msui-params-help-description {
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  line-height: var(--msui-leading-xs);
  margin: 0;

  & code {
    padding: 0 var(--msui-space-1);
    border-radius: var(--msui-radius-sm);
    background: var(--msui-control);
    font-family: var(--msui-font-mono);
    line-height: 1;
  }
}

.msui-params-help-help > * {
  grid-column: 1 / -1;
}

.msui-params-help-help > .msui-params-help-table {
  grid-template-columns: subgrid;
  gap: var(--msui-space-1) normal;
  display: grid;

  & > :nth-child(2n) {
    padding-inline-start: var(--msui-control-inset);
  }
}
}
`);
