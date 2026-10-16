/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from representation-quick.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'head': 'msp-viewer-representation-quick-head',
  'name': 'msp-viewer-representation-quick-name',
  'one': 'msp-viewer-representation-quick-one',
  'quick': 'msp-viewer-representation-quick-quick',
  'ranges': 'msp-viewer-representation-quick-ranges',
  'segment': 'msp-viewer-representation-quick-segment',
  'segments': 'msp-viewer-representation-quick-segments',
  'setting': 'msp-viewer-representation-quick-setting',
  'tile': 'msp-viewer-representation-quick-tile',
  'tiles': 'msp-viewer-representation-quick-tiles',
  'unpainted': 'msp-viewer-representation-quick-unpainted',
  'value': 'msp-viewer-representation-quick-value',
} as const;

defineStyles('viewer/views/structure/representation-quick.module.css', `@layer msp.components {
.msp-viewer-representation-quick-quick {
  gap: var(--msp-space-3);
  margin-bottom: var(--msp-space-2);
  font-size: var(--msp-font-size-sm);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-representation-quick-setting {
  gap: var(--msp-space-2);
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msp-viewer-representation-quick-ranges {
  gap: var(--msp-space-1);
  grid-template-columns: minmax(0, 1fr);
  display: grid;

  &:empty {
    display: none;
  }
}

.msp-viewer-representation-quick-head {
  align-items: baseline;
  gap: var(--msp-space-2);
  font-size: var(--msp-font-size-xs);
  display: flex;
}

.msp-viewer-representation-quick-name {
  color: var(--msp-color-text-subtle);
}

.msp-viewer-representation-quick-value {
  min-width: 0;
  color: var(--msp-color-text);
  white-space: nowrap;
  font-weight: 500;
  overflow: hidden;
}

.msp-viewer-representation-quick-segments {
  border-radius: var(--msp-radius-md);
  background: var(--msp-color-control);
  padding: 2px;
  display: flex;
}

.msp-viewer-representation-quick-segment {
  border-radius: var(--msp-radius-sm);
  min-width: 0;
  height: 30px;
  color: var(--msp-color-text-muted);
  flex: 1 1 0;
  place-items: center;
  display: grid;

  &:hover {
    color: var(--msp-color-text);
  }

  &[aria-checked="true"] {
    background: var(--msp-color-raised);
    box-shadow: 0 0 0 1px var(--msp-glass-edge),
      inset 0 0 0 1px var(--msp-selected-edge);
    color: var(--msp-color-accent);
  }
}

.msp-viewer-representation-quick-tiles {
  gap: var(--msp-space-1);
  display: flex;
}

.msp-viewer-representation-quick-tile {
  border-radius: var(--msp-radius-sm);
  background: var(--msp-color-control);
  min-width: 0;
  max-width: 44px;
  height: 22px;
  box-shadow: inset 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 88%);
  transition: transform var(--msp-motion-fast) var(--msp-ease);
  flex: 1 1 0;
  position: relative;

  &:hover {
    transform: translateY(-1px);
  }

  &[aria-checked="true"] {
    box-shadow: 0 0 0 1.5px var(--msp-color-bg),
      0 0 0 3px var(--msp-color-text);
  }
}

.msp-viewer-representation-quick-unpainted {
  background: repeating-linear-gradient(135deg,
    var(--msp-color-control) 0 4px,
    color-mix(in oklab, var(--msp-color-text-muted), transparent 75%) 4px 8px);
}

.msp-viewer-representation-quick-one {
  background: var(--msp-color-bg);
  width: 6px;
  height: 6px;
  box-shadow: 0 0 0 1px color-mix(in oklab, var(--msp-color-text), transparent 70%);
  border-radius: 50%;
  margin: -3px 0 0 -3px;
  position: absolute;
  top: 50%;
  left: 50%;
}
}
`);
