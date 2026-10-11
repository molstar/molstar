/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sequence-strip.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'chevron': 'msp-sequence-sequence-strip-chevron',
  'choice': 'msp-sequence-sequence-strip-choice',
  'gap': 'msp-sequence-sequence-strip-gap',
  'header': 'msp-sequence-sequence-strip-header',
  'info': 'msp-sequence-sequence-strip-info',
  'popover': 'msp-sequence-sequence-strip-popover',
  'residue': 'msp-sequence-sequence-strip-residue',
  'residues': 'msp-sequence-sequence-strip-residues',
  'selected': 'msp-sequence-sequence-strip-selected',
  'selector': 'msp-sequence-sequence-strip-selector',
  'spacer': 'msp-sequence-sequence-strip-spacer',
  'strip': 'msp-sequence-sequence-strip-strip',
  'value': 'msp-sequence-sequence-strip-value',
} as const;

defineStyles('views/sequence/sequence-strip.module.css', `@layer msp.components {
.msp-sequence-sequence-strip-strip {
  flex-direction: column;
  min-width: 0;
  height: 100%;
  display: flex;
}

.msp-sequence-sequence-strip-header {
  align-items: center;
  gap: var(--msp-space-1);
  height: 26px;
  padding: 0 var(--msp-space-1) 0 var(--msp-space-2);
  display: flex;
}

.msp-sequence-sequence-strip-selector {
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

.msp-sequence-sequence-strip-choice {
  min-width: 0;
  display: inline-flex;
  position: relative;
}

.msp-sequence-sequence-strip-value {
  text-overflow: ellipsis;
  white-space: nowrap;
  min-width: 0;
  max-width: 220px;
  overflow: hidden;
}

.msp-sequence-sequence-strip-chevron {
  color: var(--msp-color-text-subtle);
  flex: none;
}

.msp-sequence-sequence-strip-popover {
  z-index: var(--msp-z-overlay);
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-panel);
  min-width: 180px;
  max-width: 320px;
  box-shadow: var(--msp-glass-shadow),
    0 6px 16px #0000001f;
  position: fixed;
}

.msp-sequence-sequence-strip-info {
  color: var(--msp-color-text-muted);
  font-weight: 400;
}

.msp-sequence-sequence-strip-spacer {
  flex: 1;
}

.msp-sequence-sequence-strip-residues {
  min-width: 0;
  min-height: 0;
  padding: var(--msp-space-1) var(--msp-space-3);
  font-family: var(--msp-font-mono);
  flex-wrap: wrap;
  flex: 1;
  align-content: flex-start;
  font-size: 13px;
  line-height: 1;
  display: flex;
  overflow-y: auto;
}

.msp-sequence-sequence-strip-residue {
  text-align: center;
  cursor: pointer;
  user-select: none;
  width: 1ch;
  height: 28px;
  padding: 13px 0 2px;
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

.msp-sequence-sequence-strip-selected {
  background: color-mix(in oklab, var(--msp-color-accent), transparent 75%);
  color: var(--msp-color-accent-strong);
}

.msp-sequence-sequence-strip-gap {
  width: 1ch;
}
}
`);
