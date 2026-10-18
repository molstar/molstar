/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from story-tab.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '@molstar/prism-ui/styles/registry';

export const classes = {
  'caption': 'msp-stories-story-tab-caption',
  'captionContent': 'msp-stories-story-tab-captionContent',
  'captionScroll': 'msp-stories-story-tab-captionScroll',
  'empty': 'msp-stories-story-tab-empty',
  'header': 'msp-stories-story-tab-header',
  'number': 'msp-stories-story-tab-number',
  'scene': 'msp-stories-story-tab-scene',
  'sceneTitle': 'msp-stories-story-tab-sceneTitle',
  'scenes': 'msp-stories-story-tab-scenes',
  'tab': 'msp-stories-story-tab-tab',
  'title': 'msp-stories-story-tab-title',
} as const;

defineStyles('stories/views/story-tab.module.css', `@layer msp.components {
.msp-stories-story-tab-tab {
  gap: var(--msp-space-2);
  height: 100%;
  min-height: 0;
  padding: 0 var(--msp-panel-inset) var(--msp-panel-inset);
  font-size: var(--msp-font-size-sm);
  flex-direction: column;
  flex: 1 1 0;
  display: flex;
}

.msp-stories-story-tab-header {
  align-items: center;
  gap: var(--msp-space-2);
  min-height: 28px;
  padding-left: var(--msp-space-3);
  display: flex;
}

.msp-stories-story-tab-title {
  white-space: nowrap;
  flex: 1;
  min-width: 0;
  font-weight: 600;
  overflow: hidden;
}

.msp-stories-story-tab-scenes {
  flex: 1 1 0;
  min-height: 80px;
}

.msp-stories-story-tab-scene {
  align-items: center;
  gap: var(--msp-space-3);
  height: 28px;
  padding: 0 var(--msp-inset);
  border-radius: var(--msp-radius-md);
  display: flex;

  &:hover {
    background: var(--msp-color-hover);
  }

  &[aria-selected="true"] {
    background: var(--msp-color-accent-muted);
  }
}

.msp-stories-story-tab-number {
  min-width: 16px;
  color: var(--msp-color-text-subtle);
  font-variant-numeric: tabular-nums;
  flex: none;
}

.msp-stories-story-tab-sceneTitle {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msp-fade)), transparent);
  flex: 1;
  animation-name: msp-fade-overflow;
  animation-timing-function: linear;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msp-stories-story-tab-caption {
  border-top: 1px solid var(--msp-glass-edge);
  flex: none;
  max-height: 40%;
  display: grid;
}

.msp-stories-story-tab-captionScroll {
  min-height: 0;
}

.msp-stories-story-tab-captionContent {
  padding: var(--msp-space-3);
}

.msp-stories-story-tab-empty {
  color: var(--msp-color-text-subtle);
  margin: 0;
}
}
`);
