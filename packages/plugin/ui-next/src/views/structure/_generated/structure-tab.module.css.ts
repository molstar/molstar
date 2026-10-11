/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from structure-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'actions': 'msp-structure-structure-tab-actions',
  'chevron': 'msp-structure-structure-tab-chevron',
  'fact': 'msp-structure-structure-tab-fact',
  'facts': 'msp-structure-structure-tab-facts',
  'inspection': 'msp-structure-structure-tab-inspection',
  'inspector': 'msp-structure-structure-tab-inspector',
  'inspectorEmpty': 'msp-structure-structure-tab-inspectorEmpty',
  'inspectorHeader': 'msp-structure-structure-tab-inspectorHeader',
  'inspectorKind': 'msp-structure-structure-tab-inspectorKind',
  'inspectorTitle': 'msp-structure-structure-tab-inspectorTitle',
  'preset': 'msp-structure-structure-tab-preset',
  'presetCaption': 'msp-structure-structure-tab-presetCaption',
  'presetCheck': 'msp-structure-structure-tab-presetCheck',
  'presetGroup': 'msp-structure-structure-tab-presetGroup',
  'presetOptions': 'msp-structure-structure-tab-presetOptions',
  'presetOptionsToggle': 'msp-structure-structure-tab-presetOptionsToggle',
  'presets': 'msp-structure-structure-tab-presets',
  'quickStyles': 'msp-structure-structure-tab-quickStyles',
  'tab': 'msp-structure-structure-tab-tab',
  'tree': 'msp-structure-structure-tab-tree',
} as const;

defineStyles('views/structure/structure-tab.module.css', `@layer msp.components {
.msp-structure-structure-tab-tab {
  gap: var(--msp-space-2);
  height: 100%;
  padding: 0 var(--msp-space-2) var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  display: flex;
}

.msp-structure-structure-tab-tree {
  flex: 1 1 0;
  min-height: 80px;
  overflow: auto;
}

.msp-structure-structure-tab-actions {
  gap: var(--msp-space-1);
  flex: none;
  display: flex;
}

.msp-structure-structure-tab-inspector {
  height: clamp(180px, 42%, 420px);
  margin: 0 calc(-1 * var(--msp-space-2));
  border-top: 1px solid var(--msp-glass-edge);
  flex: none;
  overflow: auto;
}

.msp-structure-structure-tab-inspection {
  align-content: start;
  gap: var(--msp-space-2);
  padding: var(--msp-space-3) var(--msp-space-4) var(--msp-space-3);
  display: grid;
}

.msp-structure-structure-tab-inspectorEmpty {
  height: 100%;
  padding: var(--msp-space-4);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-sm);
  text-align: center;
  place-items: center;
  margin: 0;
  display: grid;
}

.msp-structure-structure-tab-inspectorHeader {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 26px;
  display: flex;
}

.msp-structure-structure-tab-inspectorTitle {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  min-width: 0;
  font-weight: 600;
  overflow: hidden;
}

.msp-structure-structure-tab-inspectorKind {
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
  text-transform: capitalize;
}

.msp-structure-structure-tab-facts {
  gap: var(--msp-space-1) var(--msp-space-4);
  font-size: var(--msp-font-size-xs);
  flex-wrap: wrap;
  margin: 0;
  display: flex;
}

.msp-structure-structure-tab-fact {
  gap: var(--msp-space-2);
  display: flex;

  & dt {
    color: var(--msp-color-text-subtle);
  }

  & dd {
    font-variant-numeric: tabular-nums;
    margin: 0;
  }
}

.msp-structure-structure-tab-quickStyles {
  gap: var(--msp-space-2);
  padding: var(--msp-space-2) var(--msp-space-2) var(--msp-space-3);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 50%);
  flex: none;
  display: grid;
}

.msp-structure-structure-tab-presets {
  gap: var(--msp-space-2);
  max-height: 50%;
  padding: var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  background: color-mix(in oklab, var(--msp-color-control), transparent 50%);
  flex: none;
  display: grid;
  overflow: auto;
}

.msp-structure-structure-tab-presetGroup {
  display: grid;
}

.msp-structure-structure-tab-presetCaption {
  padding: var(--msp-space-1) var(--msp-space-2);
  color: var(--msp-color-text-subtle);
  font-size: var(--msp-font-size-xs);
}

.msp-structure-structure-tab-preset {
  height: 28px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text);
  text-align: left;
  justify-content: space-between;
  align-items: center;
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-pressed="true"] {
    color: var(--msp-color-accent);
    font-weight: 600;
  }
}

.msp-structure-structure-tab-presetCheck {
  visibility: hidden;

  [aria-pressed="true"] > & {
    visibility: visible;
  }
}

.msp-structure-structure-tab-presetOptions {
  gap: var(--msp-space-1);
  padding-top: var(--msp-space-1);
  border-top: 1px solid var(--msp-glass-edge);
  display: grid;
}

.msp-structure-structure-tab-presetOptionsToggle {
  align-items: center;
  gap: var(--msp-space-2);
  height: 28px;
  padding: 0 var(--msp-space-2);
  border-radius: var(--msp-radius-md);
  color: var(--msp-color-text-muted);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }
}

.msp-structure-structure-tab-chevron {
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  margin-left: auto;

  [aria-expanded="false"] > & {
    transform: rotate(-90deg);
  }
}
}
`);
