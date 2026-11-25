/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Generated from playback.module.css by scripts/css.mjs; do not edit.

import { defineStyles } from '../../../styles/registry.js';

export const classes = {
  'bar': 'msui-playback-bar',
  'bare': 'msui-playback-bare',
  'card': 'msui-playback-card',
  'counter': 'msui-playback-counter',
  'deck': 'msui-playback-deck',
  'description': 'msui-playback-description',
  'descriptionContent': 'msui-playback-descriptionContent',
  'dot': 'msui-playback-dot',
  'dots': 'msui-playback-dots',
  'label': 'msui-playback-label',
  'lead': 'msui-playback-lead',
  'list': 'msui-playback-list',
  'listIndex': 'msui-playback-listIndex',
  'listItems': 'msui-playback-listItems',
  'listRow': 'msui-playback-listRow',
  'listRows': 'msui-playback-listRows',
  'listRowsContent': 'msui-playback-listRowsContent',
  'listTitle': 'msui-playback-listTitle',
  'plain': 'msui-playback-plain',
  'position': 'msui-playback-position',
  'rowWrap': 'msui-playback-rowWrap',
  'scrubber': 'msui-playback-scrubber',
  'speed': 'msui-playback-speed',
  'status': 'msui-playback-status',
  'title': 'msui-playback-title',
} as const;

defineStyles('components/playback/playback.module.css', `@layer msui.components {
.msui-playback-deck, .msui-playback-card {
  --msui-icon-size: var(--msui-icon-sm);
  border-radius: var(--msui-radius-lg);
  background: var(--msui-surface-raised);
  box-shadow: inset 0 0 0 1px var(--msui-line),
    var(--msui-shadow-float);
  color: var(--msui-text);
  font: var(--msui-text-sm) / var(--msui-leading-sm) var(--msui-font-sans);
}

.msui-playback-deck {
  grid-template-columns: minmax(0, 1fr);
  justify-self: start;
  width: fit-content;
  max-width: 100%;
  display: grid;
  overflow: clip;

  &:not(:has([data-row])) {
    display: none;
  }

  & [data-row] {
    box-shadow: 0 -1px 0 var(--msui-line);
  }
}

.msui-playback-card {
  flex-direction: column;
  min-width: 0;
  min-height: 0;
  display: flex;
  overflow: clip;
}

.msui-playback-bare {
  --msui-icon-size: var(--msui-icon-sm);
  flex-direction: column;
  flex: 1 1 0;
  min-width: 0;
  min-height: 0;
  display: flex;

  & .msui-playback-description {
    flex: 1 1 0;
    max-height: none;
  }
}

.msui-playback-rowWrap {
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msui-playback-bar {
  min-width: 0;
  height: var(--msui-header);
  padding-inline: var(--msui-surface-pad);
  align-items: center;
  gap: 0;
  display: flex;

  & button {
    color: var(--msui-text-muted);

    &:hover:not(:disabled) {
      color: var(--msui-text);
    }
  }
}

.msui-playback-lead {
  width: var(--msui-control-height);
  height: var(--msui-control-height);
  color: var(--msui-text-muted);
  flex: none;
  justify-content: center;
  align-items: center;
  display: inline-flex;
}

.msui-playback-label, .msui-playback-title {
  min-width: 0;
  max-width: 168px;
  padding-inline: var(--msui-space-2);
  white-space: nowrap;
  font-weight: var(--msui-weight-strong);
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  flex: 0 auto;
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-playback-title {
  max-width: none;
  flex: auto;
  padding-inline-start: var(--msui-inner);
}

.msui-playback-counter, .msui-playback-status {
  min-width: 40px;
  padding-inline: var(--msui-space-1);
  color: var(--msui-text-muted);
  font-size: var(--msui-text-xs);
  font-variant-numeric: tabular-nums;
  text-align: center;
  white-space: nowrap;
  flex: none;
}

.msui-playback-status {
  text-align: start;
  min-width: 56px;
  margin-inline-end: auto;

  &:empty {
    min-width: 0;
    padding: 0;
  }
}

.msui-playback-scrubber {
  min-width: 40px;
  height: var(--msui-icon);
  margin: 0 var(--msui-space-2);
  -webkit-appearance: none;
  appearance: none;
  cursor: pointer;
  background: none;
  flex: 1 0 96px;

  &::-webkit-slider-runnable-track {
    background: linear-gradient(to right,
      var(--msui-text-muted) var(--msui-progress, 0%),
      var(--msui-line-strong) var(--msui-progress, 0%));
    border-radius: 1px;
    height: 2px;
  }

  &::-webkit-slider-thumb {
    background: var(--msui-text);
    -webkit-appearance: none;
    appearance: none;
    border-radius: 50%;
    width: 10px;
    height: 10px;
    margin-top: -4px;
  }

  &::-moz-range-track {
    background: var(--msui-line-strong);
    border-radius: 1px;
    height: 2px;
  }

  &::-moz-range-progress {
    background: var(--msui-text-muted);
    border-radius: 1px;
    height: 2px;
  }

  &::-moz-range-thumb {
    background: var(--msui-text);
    border: 0;
    border-radius: 50%;
    width: 10px;
    height: 10px;
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: 2px;
    border-radius: var(--msui-radius-sm);
  }

  &:disabled {
    cursor: default;
    opacity: .5;
  }
}

.msui-playback-speed {
  flex: none;
  margin-inline-start: var(--msui-space-1);
  display: inline-flex;
}

.msui-playback-dots {
  margin-inline: var(--msui-space-1) auto;
  align-items: center;
  display: flex;
}

.msui-playback-dot {
  width: var(--msui-space-3);
  height: var(--msui-control-height);
  cursor: pointer;
  background: none;
  border: 0;
  place-items: center;
  padding: 0;
  display: grid;

  &:before {
    content: "";
    background: var(--msui-line-strong);
    border-radius: 50%;
    width: 6px;
    height: 6px;
  }

  &:hover:before {
    background: var(--msui-text-muted);
  }

  &[aria-current]:before {
    background: var(--msui-text);
  }

  &:focus-visible {
    outline: 2px solid var(--msui-focus);
    outline-offset: -2px;
    border-radius: var(--msui-radius-sm);
  }
}

.msui-playback-position {
  height: var(--msui-control-height);
  padding: 0 var(--msui-space-2);
  border-radius: var(--msui-radius-md);
  font-size: var(--msui-text-xs);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
  cursor: pointer;
  background: none;
  border: 0;
  flex: none;
  margin-inline-end: auto;

  &:hover, &[aria-expanded="true"] {
    background: var(--msui-control-hover);
  }
}

.msui-playback-title ~ .msui-playback-position, .msui-playback-title ~ .msui-playback-dots {
  margin-inline-end: 0;
}

.msui-playback-list {
  gap: var(--msui-space-1);
  min-width: 240px;
  padding: 0 var(--msui-surface-pad) var(--msui-surface-pad);
  display: grid;
}

.msui-playback-listRows {
  max-height: 216px;
}

.msui-playback-listRowsContent {
  grid-template-columns: minmax(0, 1fr);
  display: grid;
}

.msui-playback-listItems {
  grid-template-columns: minmax(0, 1fr);
  margin: 0;
  padding: 0;
  list-style: none;
  display: grid;
}

.msui-playback-listRow {
  grid-template-columns: var(--msui-icon-slot) minmax(0, 1fr);
  align-items: center;
  column-gap: var(--msui-icon-gap);
  width: 100%;
  height: var(--msui-row);
  padding: 0 var(--msui-inner);
  border-radius: var(--msui-radius-md);
  color: var(--msui-text-muted);
  font: inherit;
  text-align: start;
  cursor: pointer;
  background: none;
  border: 0;
  display: grid;

  &:hover {
    background: var(--msui-control-hover);
    color: var(--msui-text);
  }

  &[aria-current] {
    color: var(--msui-accent);
  }
}

.msui-playback-listIndex {
  color: var(--msui-text-subtle);
  font-size: var(--msui-text-xs);
  font-variant-numeric: tabular-nums;
  text-align: end;
}

.msui-playback-listTitle {
  white-space: nowrap;
  min-width: 0;
  mask-image: linear-gradient(to right, #000 calc(100% - var(--msui-fade)), transparent);
  animation-name: msui-fade-overflow;
  animation-fill-mode: both;
  animation-timeline: scroll(self inline);
  overflow: hidden;
}

.msui-playback-description {
  max-height: 40cqh;
}

.msui-playback-descriptionContent {
  min-width: 0;
  padding: 0 var(--msui-inset) var(--msui-inset);
  overflow-wrap: anywhere;
}

.msui-playback-plain {
  white-space: pre-wrap;
  margin: 0;
}
}
`);
