/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from status-line.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../styles/registry.js';

export const classes = {
  'announce': 'msui-status-status-line-announce',
  'button': 'msui-status-status-line-button',
  'chip': 'msui-status-status-line-chip',
  'chipMain': 'msui-status-status-line-chipMain',
  'chips': 'msui-status-status-line-chips',
  'countdown': 'msui-status-status-line-countdown',
  'fill': 'msui-status-status-line-fill',
  'main': 'msui-status-status-line-main',
  'pill': 'msui-status-status-line-pill',
  'pill-in': 'msui-status-status-line-pill-in',
  'progress': 'msui-status-status-line-progress',
  'spin': 'msui-status-status-line-spin',
  'status': 'msui-status-status-line-status',
  'text': 'msui-status-status-line-text',
} as const;

defineStyles('status/status-line.module.css', `@layer msui.components {
.msui-status-status-line-status {
  min-width: 0;
  max-width: 100%;
  color: var(--msui-status-ink, var(--msui-text));
  pointer-events: none;
  flex-direction: column;
  justify-content: flex-end;
  align-items: flex-start;
  display: flex;

  & > * {
    pointer-events: auto;
  }
}

.msui-status-status-line-announce {
  clip-path: inset(50%);
  white-space: nowrap;
  width: 1px;
  height: 1px;
  position: absolute;
  overflow: hidden;
}

.msui-status-status-line-pill {
  --msui-icon-size: var(--msui-icon);
  align-items: center;
  gap: var(--msui-space-1);
  max-width: 100%;
  min-height: var(--msui-row);
  font-size: var(--msui-text-sm);
  line-height: var(--msui-leading-sm);
  white-space: nowrap;
  text-shadow: 0 0 2px var(--msui-status-halo, var(--msui-bg)),
    0 0 6px var(--msui-status-halo, var(--msui-bg));
  animation: msui-status-status-line-pill-in var(--msui-enter);
  display: inline-flex;

  &[data-kind="error"] {
    color: var(--msui-danger);
  }

  &[data-kind="warning"] {
    color: var(--msui-warning);
  }

  &[data-kind="success"] {
    color: var(--msui-success);
  }

  &[data-task] {
    animation: none;
  }
}

@keyframes msui-status-status-line-pill-in {
  from {
    opacity: 0;
    transform: translateY(var(--msui-rise));
  }
}

.msui-status-status-line-main {
  align-items: center;
  gap: var(--msui-space-1);
  min-width: 0;
  height: var(--msui-row);
  padding: 0 var(--msui-space-1);
  border-radius: var(--msui-radius-md);
  color: inherit;
  font: inherit;
  text-align: start;
  background: none;
  border: 0;
  display: inline-flex;

  &:is(button) {
    cursor: pointer;
  }

  &:is(button):hover .msui-status-status-line-text {
    text-underline-offset: 3px;
    text-decoration: underline;
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}

.msui-status-status-line-text {
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-status-status-line-button {
  --msui-icon-size: var(--msui-icon-sm);
  width: var(--msui-row);
  height: var(--msui-row);
  border-radius: var(--msui-radius-md);
  color: inherit;
  opacity: .8;
  cursor: pointer;
  background: none;
  border: 0;
  flex: none;
  place-items: center;
  padding: 0;
  display: grid;
  position: relative;

  &:hover {
    background: var(--msui-control-hover);
    opacity: 1;
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }

  &[hidden] {
    display: none;
  }
}

.msui-status-status-line-countdown {
  width: var(--msui-icon-lg);
  height: var(--msui-icon-lg);
  pointer-events: none;
  position: absolute;
  top: 50%;
  left: 50%;
  transform: translate(-50%, -50%) rotate(-90deg);

  & circle {
    fill: none;
    stroke: currentColor;
    stroke-width: 1.5px;
    stroke-dasharray: 100;
    stroke-dashoffset: 0;
    opacity: .6;
    animation: msui-status-status-line-countdown var(--msui-countdown, 8s) linear forwards;
  }

  .msui-status-status-line-pill:hover &, .msui-status-status-line-pill:focus-within & {
    & circle {
      animation-play-state: paused;
    }
  }
}

@keyframes msui-status-status-line-countdown {
  to {
    stroke-dashoffset: 100px;
  }
}

.msui-status-status-line-progress {
  width: calc(12 * var(--msui-unit));
  height: var(--msui-unit);
  border-radius: var(--msui-radius-sm);
  background: var(--msui-line-strong);
  flex: none;
  display: block;
  overflow: hidden;

  &[hidden] {
    display: none;
  }
}

.msui-status-status-line-fill {
  background: var(--msui-accent);
  height: 100%;
  display: block;
}

.msui-status-status-line-spin {
  animation: 1s linear infinite msui-status-status-line-spin;
}

@keyframes msui-status-status-line-spin {
  to {
    transform: rotate(360deg);
  }
}

@media (prefers-reduced-motion: reduce) {
  .msui-status-status-line-countdown {
    display: none;
  }

  .msui-status-status-line-spin {
    animation: none;
  }
}

.msui-status-status-line-chips {
  justify-content: flex-end;
  gap: var(--msui-space-2);
  flex-wrap: wrap;
  min-width: 0;
  display: flex;
}

.msui-status-status-line-chip {
  --msui-icon-size: var(--msui-icon-sm);
  max-width: 100%;
  height: var(--msui-header);
  padding-inline: var(--msui-space-1);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: var(--msui-shadow-float);
  color: var(--msui-text);
  white-space: nowrap;
  align-items: center;
  display: inline-flex;
}

.msui-status-status-line-chipMain {
  min-width: 0;
  height: var(--msui-row);
  padding: 0 var(--msui-space-2);
  border-radius: var(--msui-radius-md);
  color: inherit;
  font: inherit;
  font-size: var(--msui-text-md);
  line-height: var(--msui-leading-md);
  background: none;
  border: 0;
  align-items: center;
  display: inline-flex;

  &:enabled {
    cursor: pointer;
  }

  &:enabled:hover {
    background: var(--msui-control-hover);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
  }
}
}
`);
