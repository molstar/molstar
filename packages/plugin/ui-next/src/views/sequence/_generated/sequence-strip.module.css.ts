/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from sequence-strip.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'chain': 'msp-sequence-sequence-strip-chain',
  'chainInfo': 'msp-sequence-sequence-strip-chainInfo',
  'chainName': 'msp-sequence-sequence-strip-chainName',
  'gap': 'msp-sequence-sequence-strip-gap',
  'residue': 'msp-sequence-sequence-strip-residue',
  'residues': 'msp-sequence-sequence-strip-residues',
  'selected': 'msp-sequence-sequence-strip-selected',
  'strip': 'msp-sequence-sequence-strip-strip',
} as const;

defineStyles('views/sequence/sequence-strip.module.css', `@layer msp.components {
.msp-sequence-sequence-strip-strip {
  align-items: stretch;
  min-width: 0;
  height: 100%;
  display: flex;
}

.msp-sequence-sequence-strip-chain {
  padding: 0 var(--msp-space-4);
  border-right: 1px solid var(--msp-color-border);
  font-size: var(--msp-font-size-sm);
  white-space: nowrap;
  flex-direction: column;
  justify-content: center;
  gap: 2px;
  display: flex;
}

.msp-sequence-sequence-strip-chainName {
  font-weight: 600;
}

.msp-sequence-sequence-strip-chainInfo {
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-xs);
}

.msp-sequence-sequence-strip-residues {
  min-width: 0;
  padding: var(--msp-space-2) var(--msp-space-3);
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
  padding: 14px 0 2px;
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
