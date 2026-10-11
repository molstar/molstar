/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from dock.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'barTab': 'msp-shell-dock-barTab',
  'body': 'msp-shell-dock-body',
  'bottom': 'msp-shell-dock-bottom',
  'collapsedBar': 'msp-shell-dock-collapsedBar',
  'compact': 'msp-shell-dock-compact',
  'content': 'msp-shell-dock-content',
  'dock': 'msp-shell-dock-dock',
  'header': 'msp-shell-dock-header',
  'left': 'msp-shell-dock-left',
  'right': 'msp-shell-dock-right',
  'splitter': 'msp-shell-dock-splitter',
  'strip': 'msp-shell-dock-strip',
  'tab': 'msp-shell-dock-tab',
  'title': 'msp-shell-dock-title',
  'top': 'msp-shell-dock-top',
} as const;

defineStyles('views/shell/dock.module.css', `@layer msp.components {
.msp-shell-dock-dock {
  background: var(--msp-color-panel);
  min-width: 0;
  min-height: 0;
  display: flex;
}

.msp-shell-dock-right {
  border-left: 1px solid var(--msp-color-border);
  flex-direction: row;
  grid-area: right;
}

.msp-shell-dock-left {
  border-right: 1px solid var(--msp-color-border);
  flex-direction: row-reverse;
  grid-area: left;
}

.msp-shell-dock-top, .msp-shell-dock-bottom {
  flex-direction: row-reverse;

  & .msp-shell-dock-content {
    width: auto;
    height: var(--msp-dock-size);
    flex: 1;
  }

  & .msp-shell-dock-header, & .msp-shell-dock-splitter {
    display: none;
  }

  & .msp-shell-dock-strip {
    border-right: 1px solid var(--msp-color-border);
  }

  & .msp-shell-dock-tab {
    width: 28px;
    height: 28px;
  }
}

.msp-shell-dock-top {
  border-bottom: 1px solid var(--msp-color-border);
  grid-area: top;
}

.msp-shell-dock-bottom {
  border-top: 1px solid var(--msp-color-border);
  grid-area: bottom;
}

.msp-shell-dock-content {
  width: var(--msp-dock-size);
  flex-direction: column;
  min-width: 0;
  min-height: 0;
  display: flex;
  position: relative;
}

.msp-shell-dock-header {
  align-items: center;
  gap: var(--msp-space-3);
  height: var(--msp-row-height);
  padding: 0 var(--msp-space-2) 0 var(--msp-space-4);
  border-bottom: 1px solid var(--msp-color-border);
  font-weight: 600;
  display: flex;
}

.msp-shell-dock-title {
  text-overflow: ellipsis;
  white-space: nowrap;
  flex: 1;
  overflow: hidden;
}

.msp-shell-dock-body {
  flex: 1;
  min-height: 0;
  overflow: auto;
}

.msp-shell-dock-strip {
  gap: var(--msp-space-1);
  padding: var(--msp-space-2);
  background: var(--msp-color-bg);
  flex-direction: column;
  display: flex;
}

.msp-shell-dock-right .msp-shell-dock-strip {
  border-left: 1px solid var(--msp-color-border);
}

.msp-shell-dock-left .msp-shell-dock-strip {
  border-right: 1px solid var(--msp-color-border);
}

.msp-shell-dock-tab {
  border-radius: var(--msp-radius-md);
  width: 36px;
  height: 36px;
  color: var(--msp-color-text-muted);
  transition: background-color var(--msp-motion-fast) var(--msp-ease),
    color var(--msp-motion-fast) var(--msp-ease);
  place-items: center;
  display: grid;

  &:hover {
    background: var(--msp-color-hover);
    color: var(--msp-color-text);
  }

  &[aria-selected="true"] {
    background: var(--msp-color-control);
    color: var(--msp-color-accent);
  }
}

.msp-shell-dock-splitter {
  z-index: 1;
  cursor: col-resize;
  touch-action: none;
  width: 6px;
  position: absolute;
  top: 0;
  bottom: 0;

  &:hover, &[data-dragging] {
    background: color-mix(in oklab, var(--msp-color-accent), transparent 60%);
  }
}

.msp-shell-dock-right .msp-shell-dock-splitter {
  left: -3px;
}

.msp-shell-dock-left .msp-shell-dock-splitter {
  right: -3px;
}

.msp-shell-dock-collapsedBar {
  align-items: center;
  gap: var(--msp-space-3);
  height: 24px;
  padding: 0 var(--msp-space-3);
  flex: 1;
  display: flex;
}

.msp-shell-dock-barTab {
  align-items: center;
  gap: var(--msp-space-1);
  color: var(--msp-color-text-muted);
  font-size: var(--msp-font-size-sm);
  display: inline-flex;

  &:hover {
    color: var(--msp-color-text);
  }
}

.msp-shell-dock-compact {
  border-left: 0;
  border-top: 1px solid var(--msp-color-border);
  flex-direction: column;

  & .msp-shell-dock-content {
    width: auto;
    height: min(45cqh, 420px);
  }

  & .msp-shell-dock-strip {
    border-left: 0;
    border-top: 1px solid var(--msp-color-border);
    flex-direction: row;
    justify-content: space-around;
  }

  & .msp-shell-dock-splitter {
    display: none;
  }
}
}
`);
