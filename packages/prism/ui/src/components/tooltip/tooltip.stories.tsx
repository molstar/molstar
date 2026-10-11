/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { onCleanup } from '../../dom/scope.js';
import { story } from '../../testing/story-root.js';
import { Icon } from '../icon/icon.js';
import { installTooltips } from './tooltip.js';

const meta: Meta = { title: 'Components/Tooltip' };
export default meta;

const panel =
  'display: grid; gap: var(--msp-space-4); width: 300px; padding: var(--msp-space-4) var(--msp-space-4) 48px; background: var(--msp-color-panel); color: var(--msp-color-text); font-size: var(--msp-font-size-sm);';
const row =
  'display: flex; gap: var(--msp-space-3); align-items: center; justify-content: center; padding: var(--msp-space-6) 0;';
const muted = 'color: var(--msp-color-text-muted);';

// Removed with the story: the ref runs in the story's scope, so the uninstall function goes to its cleanup.
const tips = (el: HTMLElement) => onCleanup(installTooltips(el));

const iconButton =
  'display: grid; place-items: center; width: 28px; height: 28px; border-radius: var(--msp-radius-md); background: var(--msp-color-control); color: var(--msp-color-text);';

/** Icon buttons with `data-tip`, placed by `data-tip-side`: bottom (the default), left and top. Hover or focus one. */
export const Sides: StoryObj = {
  render: () =>
    story(() => (
      <div style={panel} ref={tips}>
        <div style={row}>
          <button type="button" style={iconButton} aria-label="Refresh" data-tip="Bottom (default)">
            <Icon name="refresh-cw" size={16} />
          </button>
          <button type="button" style={iconButton} aria-label="Add" data-tip="Left" data-tip-side="left">
            <Icon name="plus" size={16} />
          </button>
          <button type="button" style={iconButton} aria-label="Remove" data-tip="Top" data-tip-side="top">
            <Icon name="trash-2" size={16} />
          </button>
        </div>
        <span style={muted}>Shown after a short rest; moving to a neighbour switches at once.</span>
      </div>
    )),
};

/** `data-tip-overflow`: a label shows its own text as a tooltip, only while it is cut off. */
export const Overflow: StoryObj = {
  render: () =>
    story(() => (
      <div style={panel} ref={tips}>
        <span
          data-tip-overflow
          style="display: block; width: 120px; overflow: hidden; white-space: nowrap; mask-image: linear-gradient(to right, #000 calc(100% - 20px), transparent)"
        >
          A long label that does not fit
        </span>
        <span
          data-tip-overflow
          style="display: block; width: 120px; overflow: hidden; white-space: nowrap; mask-image: linear-gradient(to right, #000 calc(100% - 20px), transparent)"
        >
          Fits
        </span>
        <span style={muted}>Only the first one is cut off, so only it has a tooltip.</span>
      </div>
    )),
};
