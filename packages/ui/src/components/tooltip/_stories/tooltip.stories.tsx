/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { onCleanup } from '../../../dom/index.js';
import { story } from '../../../_stories/support.js';
import { Icon } from '../../icon/icon.js';
import { installTooltips } from '../tooltip.js';

const meta: Meta = { title: 'Components/Tooltip' };
export default meta;

const panel =
  'display: grid; gap: var(--msui-space-3); width: 300px; padding: var(--msui-space-3) var(--msui-space-3) 48px; background: var(--msui-surface); color: var(--msui-text); font-size: var(--msui-text-sm);';
const row =
  'display: flex; gap: var(--msui-space-2); align-items: center; justify-content: center; padding: var(--msui-space-5) 0;';
const muted = 'color: var(--msui-text-muted);';

// Removed with the story: the ref runs in the story's scope, so the uninstall function goes to its cleanup.
const tips = (el: HTMLElement) => onCleanup(installTooltips(el));

const iconButton =
  'display: grid; place-items: center; width: 28px; height: 28px; border-radius: var(--msui-radius-md); background: var(--msui-control); color: var(--msui-text);';

/** Icon buttons with `data-tip`, placed by `data-tip-side`: bottom (the default), left and top. Hover or focus one. */
export const Sides: StoryObj = {
  render: () =>
    story(() => (
      <div style={panel} ref={tips}>
        <div style={row}>
          <button type="button" style={iconButton} aria-label="Refresh" data-tip="Bottom (default)">
            <Icon name="refresh-cw" />
          </button>
          <button type="button" style={iconButton} aria-label="Add" data-tip="Left" data-tip-side="left">
            <Icon name="plus" />
          </button>
          <button type="button" style={iconButton} aria-label="Remove" data-tip="Top" data-tip-side="top">
            <Icon name="trash-2" />
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
