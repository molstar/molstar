/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sequence-strip.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'chevron': 'msp-viewer-sequence-strip-chevron',
  'choice': 'msp-viewer-sequence-strip-choice',
  'gap': 'msp-viewer-sequence-strip-gap',
  'header': 'msp-viewer-sequence-strip-header',
  'info': 'msp-viewer-sequence-strip-info',
  'popover': 'msp-viewer-sequence-strip-popover',
  'residue': 'msp-viewer-sequence-strip-residue',
  'residues': 'msp-viewer-sequence-strip-residues',
  'residuesContent': 'msp-viewer-sequence-strip-residuesContent',
  'selected': 'msp-viewer-sequence-strip-selected',
  'selector': 'msp-viewer-sequence-strip-selector',
  'spacer': 'msp-viewer-sequence-strip-spacer',
  'strip': 'msp-viewer-sequence-strip-strip',
  'value': 'msp-viewer-sequence-strip-value',
} as const;

defineStyles('viewer/views/sequence/sequence-strip.module.css', `@layer msp.components {
.msp-viewer-sequence-strip-strip {
  flex-direction: column;
  flex: 1 1 0;
  min-width: 0;
  height: 100%;
  min-height: 0;
  display: flex;
}

.msp-viewer-sequence-strip-header {
  align-items: center;
  gap: var(--msp-space-1);
  height: 26px;
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  display: flex;
}

.msp-viewer-sequence-strip-selector {
  align-items: center;
  gap: var(--msp-space-2);
  min-width: 0;
  height: 22px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  font-weight: 500;
  display: inline-flex;

  &:hover, &[aria-expanded="true"] {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[data-strong] {
    color: var(--msp-color-text);
    font-weight: 600;
  }
}

.msp-viewer-sequence-strip-choice {
  min-width: 0;
  display: inline-flex;
  position: relative;
}

.msp-viewer-sequence-strip-value {
  min-width: 0;
  max-width: 220px;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  white-space: nowrap;
  animation-name: msp-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-viewer-sequence-strip-chevron {
  color: var(--msp-color-text-subtle);
  flex: none;
}

.msp-viewer-sequence-strip-popover {
  z-index: var(--msp-z-overlay);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  min-width: 180px;
  max-width: 320px;
  box-shadow: var(--msp-glass-shadow),
    0 6px 16px #0000001f;
  position: fixed;
}

.msp-viewer-sequence-strip-info {
  color: var(--msp-color-text-muted);
  font-weight: 400;
}

.msp-viewer-sequence-strip-spacer {
  flex: 1;
}

.msp-viewer-sequence-strip-residues {
  flex: 1;
  min-width: 0;
  min-height: 0;
}

.msp-viewer-sequence-strip-residuesContent {
  padding: var(--msp-space-1) var(--msp-space-3);
  font-family: var(--msp-font-mono);
  flex-wrap: wrap;
  align-content: flex-start;
  font-size: 13px;
  line-height: 1;
  display: flex;
}

.msp-viewer-sequence-strip-residue {
  height: 28px;
  padding: 13px 0 var(--msp-space-1);
  text-align: center;
  cursor: pointer;
  user-select: none;
  width: 1ch;
  position: relative;

  &:hover {
    color: var(--msp-color-accent);
    background: var(--msp-color-hover);
  }

  &[data-number]:before {
    content: attr(data-number);
    font-family: var(--msp-font-sans);
    color: var(--msp-color-text-subtle);
    font-size: 9px;
    position: absolute;
    top: 2px;
    left: 0;
  }
}

.msp-viewer-sequence-strip-selected {
  background: color-mix(in oklab, var(--msp-color-accent), transparent 75%);
  color: var(--msp-color-accent-strong);
}

.msp-viewer-sequence-strip-gap {
  width: 1ch;
}
}
`);
