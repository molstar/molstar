/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from canvas.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'bottom': 'msp-shell-canvas-bottom',
  'card': 'msp-shell-canvas-card',
  'card-in': 'msp-shell-canvas-card-in',
  'cardBody': 'msp-shell-canvas-cardBody',
  'cardHeader': 'msp-shell-canvas-cardHeader',
  'cardTitle': 'msp-shell-canvas-cardTitle',
  'message': 'msp-shell-canvas-message',
  'pill': 'msp-shell-canvas-pill',
  'pill-in': 'msp-shell-canvas-pill-in',
  'pillButton': 'msp-shell-canvas-pillButton',
  'progress': 'msp-shell-canvas-progress',
  'separator': 'msp-shell-canvas-separator',
  'spin': 'msp-shell-canvas-spin',
  'status': 'msp-shell-canvas-status',
  'text': 'msp-shell-canvas-text',
  'timeline': 'msp-shell-canvas-timeline',
  'tool': 'msp-shell-canvas-tool',
  'toolbar': 'msp-shell-canvas-toolbar',
} as const;

defineStyles('views/shell/canvas.module.css', `@layer msp.components {
.msp-shell-canvas-toolbar {
  top: var(--msp-space-4);
  right: var(--msp-space-4);
  gap: var(--msp-space-1);
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  border: var(--msp-glass-border);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  flex-direction: column;
  display: flex;
  position: absolute;
}

.msp-shell-canvas-tool {
  border-radius: var(--msp-radius-md);
  width: 34px;
  height: 34px;
  color: var(--msp-color-text-muted);
  transition: background-color var(--msp-motion-fast) var(--msp-ease),
    color var(--msp-motion-fast) var(--msp-ease);
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[aria-pressed="true"] {
    background: var(--msp-color-control);
    color: var(--msp-color-accent);
  }
}

.msp-shell-canvas-separator {
  height: 1px;
  margin: var(--msp-space-1) var(--msp-space-2);
  background: var(--msp-color-border);
}

.msp-shell-canvas-card {
  right: calc(var(--msp-space-4) + 34px + 2 * var(--msp-space-2) + var(--msp-space-3));
  z-index: var(--msp-z-overlay);
  width: min(360px, 100% - 96px);
  max-height: calc(100% - 2 * var(--msp-space-4));
  border-radius: var(--msp-radius-lg);
  animation: msp-shell-canvas-card-in var(--msp-motion-base) var(--msp-ease);
  background: var(--msp-glass-bg);
  border: var(--msp-glass-border);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  flex-direction: column;
  display: flex;
  position: absolute;
  overflow: hidden;
}

@keyframes msp-shell-canvas-card-in {
  from {
    opacity: 0;
    transform: translateX(8px);
  }
}

.msp-shell-canvas-cardHeader {
  align-items: center;
  gap: var(--msp-space-3);
  height: var(--msp-row-height);
  padding: 0 var(--msp-space-2) 0 var(--msp-space-4);
  border-bottom: 1px solid var(--msp-color-border);
  font-weight: 600;
  display: flex;
}

.msp-shell-canvas-cardTitle {
  flex: 1;
}

.msp-shell-canvas-cardBody {
  min-height: 0;
  overflow: auto;
}

.msp-shell-canvas-bottom {
  left: var(--msp-space-4);
  right: var(--msp-space-4);
  bottom: var(--msp-space-4);
  align-items: end;
  gap: var(--msp-space-3);
  pointer-events: none;
  grid-template-columns: minmax(0, 1fr) auto minmax(0, 1fr);
  grid-template-areas: "status timeline .";
  display: grid;
  position: absolute;

  & > * {
    pointer-events: auto;
  }
}

@container canvas (width < 960px) {
  .msp-shell-canvas-bottom {
    grid-template-columns: minmax(0, 1fr);
    grid-template-areas: "timeline"
                         "status";
  }
}

.msp-shell-canvas-status {
  align-items: center;
  gap: var(--msp-space-2);
  pointer-events: none;
  flex-wrap: wrap;
  grid-area: status;
  display: flex;

  & > * {
    pointer-events: auto;
  }
}

.msp-shell-canvas-timeline {
  grid-area: timeline;
  justify-self: center;
}

.msp-shell-canvas-pill {
  align-items: center;
  gap: var(--msp-space-2);
  max-width: 100%;
  height: 26px;
  padding: 0 var(--msp-space-3);
  color: var(--msp-color-text);
  font-size: var(--msp-font-size-sm);
  white-space: nowrap;
  background: var(--msp-glass-bg);
  border: var(--msp-glass-border);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
  border-radius: 13px;
  display: inline-flex;
}

.msp-shell-canvas-message {
  animation: msp-shell-canvas-pill-in var(--msp-motion-base) var(--msp-ease);

  &[data-kind="error"] {
    color: var(--msp-color-danger);
  }

  &[data-kind="warning"] {
    color: var(--msp-color-warning);
  }

  &[data-kind="success"] {
    color: var(--msp-color-success);
  }
}

@keyframes msp-shell-canvas-pill-in {
  from {
    opacity: 0;
    transform: translateY(6px);
  }
}

.msp-shell-canvas-text {
  text-overflow: ellipsis;
  overflow: hidden;
}

.msp-shell-canvas-progress {
  background: var(--msp-color-border);
  border-radius: 2px;
  width: 48px;
  height: 4px;
  overflow: hidden;

  & > div {
    background: var(--msp-color-accent);
    height: 100%;
  }
}

.msp-shell-canvas-spin {
  animation: 1s linear infinite msp-shell-canvas-spin;
}

@keyframes msp-shell-canvas-spin {
  to {
    transform: rotate(360deg);
  }
}

.msp-shell-canvas-pillButton {
  margin-right: calc(-1 * var(--msp-space-2));
  color: var(--msp-color-text-muted);
  border-radius: 50%;
  place-items: center;
  display: grid;

  &:hover {
    color: var(--msp-color-text);
  }
}
}
`);
