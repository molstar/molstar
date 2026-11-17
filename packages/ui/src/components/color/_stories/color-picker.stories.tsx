/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { derive } from '@molstar/core/reactive';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { onCleanup } from '../../../dom/index.js';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { ColorPicker, colorToHex, colorToRgb } from '../color-picker.js';

const meta: Meta = { title: 'Components/Color picker' };
export default meta;

/**
 * In a form, under its row (as a select's options open): the smooth palette to click or drag, the palette's hues, the
 * weights of the color picked, and a field that takes a color as people write one.
 */
export const InAForm: StoryObj = {
  render: () =>
    story(() => {
      const color = state(0x0090ff);
      return (
        <div
          ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
          class="msui-surface"
          style="width: 320px; margin: var(--msui-space-4);"
        >
          <div class="msui-form">
            <span class="msui-label">Color</span>
            <span class="msui-value" style="display: flex; align-items: center; gap: var(--msui-icon-gap);">
              <span
                style={{
                  width: 'var(--msui-icon)',
                  height: 'var(--msui-icon)',
                  borderRadius: 'var(--msui-radius-sm)',
                  background: derive(color, colorToHex),
                }}
              />
              {derive(color, colorToRgb)}
            </span>
            <div style="grid-column: 1 / -1;">
              <ColorPicker value={color} onChange={(c) => color.set(c)} />
            </div>
          </div>
        </div>
      );
    }),
};
