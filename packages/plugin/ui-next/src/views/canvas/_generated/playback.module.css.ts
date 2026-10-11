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
  'progress': 'msp-canvas-playback-progress',
  'scrubber': 'msp-canvas-playback-scrubber',
  'title': 'msp-canvas-playback-title',
} as const;

defineStyles('views/canvas/playback.module.css', `@layer msp.components {
.msp-canvas-playback-cluster {
  top: var(--msp-space-4);
  left: var(--msp-space-4);
  gap: var(--msp-space-2);
  width: min(320px, 100% - 96px);
  display: grid;
  position: absolute;
}

.msp-canvas-playback-card {
  border-radius: var(--msp-radius-lg);
  background: color-mix(in oklab, var(--msp-color-panel), transparent 12%);
  box-shadow: var(--msp-shadow-float);
  backdrop-filter: blur(6px);
}

.msp-canvas-playback-bar {
  align-items: center;
  gap: var(--msp-space-1);
  height: 34px;
  padding: 0 var(--msp-space-2);
  display: flex;
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

.msp-canvas-playback-progress {
  height: 2px;
  margin: 0 var(--msp-space-3);
  background: var(--msp-color-accent);
  transform-origin: 0;
  transition: transform var(--msp-motion-fast) linear;
}

.msp-canvas-playback-description {
  max-height: 40cqh;
  padding: var(--msp-space-3) var(--msp-space-4) var(--msp-space-4);
  border-top: 1px solid var(--msp-color-border);
  font-size: var(--msp-font-size-sm);
  overflow: auto;
}

.msp-canvas-playback-scrubber {
  min-width: 40px;
  accent-color: var(--msp-color-accent);
  flex: 1;
}
}
`);
