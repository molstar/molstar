/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from action-footer.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'apply': 'msui-action-footer-apply',
  'footer': 'msui-action-footer-footer',
  'info': 'msui-action-footer-info',
  'reset': 'msui-action-footer-reset',
  'secondary': 'msui-action-footer-secondary',
  'spin': 'msui-action-footer-spin',
} as const;

defineStyles('components/action-footer/action-footer.module.css', `@layer msui.components {
.msui-action-footer-footer {
  height: var(--msui-row);
  padding-inline: var(--msui-inner);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
  --msui-icon-size: var(--msui-icon);
  align-items: center;
  display: flex;
}

.msui-action-footer-reset, .msui-action-footer-secondary, .msui-action-footer-apply {
  height: var(--msui-control-height);
  border-radius: var(--msui-radius-md);
  font: inherit;
  cursor: pointer;
  background: none;
  border: 0;
  align-items: center;
  display: flex;

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &:disabled {
    color: var(--msui-text-subtle);
    cursor: default;
  }
}

.msui-action-footer-reset {
  width: var(--msui-control-height);
  margin-inline: calc((var(--msui-icon-slot) - var(--msui-control-height)) / 2);
  color: var(--msui-text-subtle);
  flex: none;
  justify-content: center;
  padding: 0;

  &:hover, &:focus-visible {
    color: var(--msui-text);
  }
}

.msui-action-footer-info {
  min-width: 0;
  color: var(--msui-text-subtle);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  flex: 1;
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-action-footer-reset + .msui-action-footer-info {
  margin-inline-start: var(--msui-icon-gap);
}

.msui-action-footer-apply {
  gap: var(--msui-icon-gap);
  padding-inline: var(--msui-control-inset) var(--msui-space-1);
  color: var(--msui-accent);
  font-weight: var(--msui-weight-strong);
  flex: none;
  margin-inline-start: auto;

  &:hover:not(:disabled) {
    color: var(--msui-accent-strong);

    & > span {
      text-underline-offset: 3px;
      text-decoration: underline;
    }
  }
}

.msui-action-footer-secondary {
  padding-inline: var(--msui-control-inset);
  color: var(--msui-text-muted);
  flex: none;

  &:hover:not(:disabled) {
    color: var(--msui-text);

    & > span {
      text-underline-offset: 3px;
      text-decoration: underline;
    }
  }
}

.msui-action-footer-footer:not(:has(.msui-action-footer-info)) .msui-action-footer-secondary {
  margin-inline-start: auto;
}

.msui-action-footer-secondary ~ .msui-action-footer-apply {
  margin-inline-start: 0;
}

.msui-action-footer-spin {
  animation: .9s linear infinite msui-action-footer-spin;
}

@keyframes msui-action-footer-spin {
  to {
    transform: rotate(360deg);
  }
}

@media (prefers-reduced-motion: reduce) {
  .msui-action-footer-spin {
    animation-duration: 2.7s;
  }
}
}
`);
