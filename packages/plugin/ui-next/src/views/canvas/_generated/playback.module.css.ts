/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from playback.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'bar': 'msp-canvas-playback-bar',
  'card': 'msp-canvas-playback-card',
  'cluster': 'msp-canvas-playback-cluster',
  'counter': 'msp-canvas-playback-counter',
  'description': 'msp-canvas-playback-description',
  'scrubber': 'msp-canvas-playback-scrubber',
  'timeline': 'msp-canvas-playback-timeline',
  'title': 'msp-canvas-playback-title',
} as const;

defineStyles('views/canvas/playback.module.css', `@layer msp.components {
.msp-canvas-playback-cluster {
  gap: var(--msp-space-2);
  display: grid;
}

.msp-canvas-playback-card {
  border-radius: var(--msp-radius-lg);
  background: var(--msp-glass-bg);
  box-shadow: var(--msp-glass-shadow);
  -webkit-backdrop-filter: var(--msp-glass-filter);
}

.msp-canvas-playback-bar {
  align-items: center;
  gap: var(--msp-space-1);
  height: 34px;
  padding: 0 var(--msp-space-2);
  display: flex;

  & button {
    color: var(--msp-color-text-muted);

    &:hover:not(:disabled) {
      color: var(--msp-color-text);
    }

    &[aria-pressed="true"] {
      color: var(--msp-color-accent);
    }
  }

  & svg {
    width: 16px;
    height: 16px;
  }
}

.msp-canvas-playback-counter {
  min-width: 44px;
  font-size: var(--msp-font-size-sm);
  font-variant-numeric: tabular-nums;
  text-align: center;
  color: var(--msp-color-text-muted);
}

.msp-canvas-playback-title {
  min-width: 0;
  padding: 0 var(--msp-space-2);
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  font-weight: 600;
  overflow: hidden;
}

.msp-canvas-playback-description {
  max-height: 40cqh;
  padding: var(--msp-space-3) var(--msp-space-4) var(--msp-space-4);
  font-size: var(--msp-font-size-sm);
  overflow: auto;
}

.msp-canvas-playback-timeline {
  width: 100%;
}

.msp-canvas-playback-scrubber {
  min-width: 40px;
  accent-color: var(--msp-color-accent);
  flex: 1;
}
}
`);
