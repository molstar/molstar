/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from audio-bar.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'action': 'msp-viewer-audio-bar-action',
  'bar': 'msp-viewer-audio-bar-bar',
  'fill': 'msp-viewer-audio-bar-fill',
  'knob': 'msp-viewer-audio-bar-knob',
  'knobTrack': 'msp-viewer-audio-bar-knobTrack',
  'play': 'msp-viewer-audio-bar-play',
  'time': 'msp-viewer-audio-bar-time',
  'track': 'msp-viewer-audio-bar-track',
} as const;

defineStyles('viewer/views/state/audio-bar.module.css', `@layer msp.components {
.msp-viewer-audio-bar-bar {
  align-items: center;
  gap: var(--msp-space-2);
  height: 32px;
  padding: 0 var(--msp-space-1);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-well);
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-viewer-audio-bar-play {
  background: var(--msp-color-accent);
  width: 24px;
  height: 24px;
  color: var(--msp-color-on-accent);
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  border-radius: 50%;
  flex: none;
  place-items: center;
  display: grid;

  &:hover {
    transform: scale(1.06);
  }
}

.msp-viewer-audio-bar-time {
  min-width: 4ch;
  color: var(--msp-color-text-subtle);
  font-variant-numeric: tabular-nums;
  text-align: center;
  flex: none;
}

.msp-viewer-audio-bar-track {
  cursor: pointer;
  touch-action: none;
  outline: none;
  flex: 1;
  min-width: 0;
  height: 16px;
  position: relative;

  &:before {
    content: "";
    background: var(--msp-color-pressed);
    border-radius: 1.5px;
    height: 3px;
    margin-top: -1.5px;
    position: absolute;
    top: 50%;
    left: 0;
    right: 0;
  }

  &:focus-visible {
    box-shadow: 0 0 0 1px var(--msp-color-focus-ring);
    border-radius: var(--msp-radius-sm);
  }
}

.msp-viewer-audio-bar-fill {
  background: var(--msp-color-accent);
  transform-origin: 0;
  will-change: transform;
  border-radius: 1.5px;
  width: 100%;
  height: 3px;
  margin-top: -1.5px;
  position: absolute;
  top: 50%;
  left: 0;
}

.msp-viewer-audio-bar-knobTrack {
  pointer-events: none;
  will-change: transform;
  position: absolute;
  inset: 0;
}

.msp-viewer-audio-bar-knob {
  background: var(--msp-color-accent);
  width: 10px;
  height: 10px;
  box-shadow: 0 0 0 2px var(--msp-color-well);
  border-radius: 50%;
  margin-top: -5px;
  position: absolute;
  top: 50%;
  left: -5px;
}

.msp-viewer-audio-bar-action {
  border-radius: var(--msp-radius-sm);
  width: 22px;
  height: 22px;
  color: var(--msp-color-text-muted);
  flex: none;
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}
}
`);
