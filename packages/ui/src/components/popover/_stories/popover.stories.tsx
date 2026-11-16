/** @jsxImportSource @molstar/ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { onCleanup, Show } from '../../../dom/index.js';
import { story } from '../../../_stories/support.js';
import { state } from '../../../state.js';
import { Button } from '../../button/button.js';
import { Icon } from '../../icon/icon.js';
import { Select } from '../../select/select.js';
import { TextInput } from '../../text-input/text-input.js';
import { installTooltips } from '../../tooltip/tooltip.js';
import { Popover } from '../popover.js';

const meta: Meta = { title: 'Components/Popover' };
export default meta;

/**
 * A header's settings: the button opens a panel of controls under it, its right edge on the button's. Its title is a
 * header on the inset; the form in it puts its labels on the inset and its fields up to the inner padding. Escape or a
 * press outside closes it; a select opened from it keeps it.
 */
export const Settings: StoryObj = {
  render: () =>
    story(() => {
      const open = state(false);
      const quality = state('auto');
      const name = state('Scene');
      let anchor!: HTMLElement;
      return (
        <div
          ref={(el: HTMLElement) => onCleanup(installTooltips(el))}
          class="msui-surface"
          style="width: 420px; margin: var(--msui-space-4);"
        >
          <div class="msui-header" style="justify-content: space-between; padding-inline-end: var(--msui-space-1);">
            Structure
            <span ref={(el: HTMLElement) => (anchor = el)}>
              <Button
                variant="ghost"
                icon={<Icon name="settings" />}
                title="Scene defaults"
                pressed={open}
                onClick={() => open.set(!open.value)}
              />
            </span>
          </div>
          <Show when={open}>
            {() => (
              <Popover anchor={anchor} title="Scene defaults" onClose={() => open.set(false)}>
                <div class="msui-form">
                  <span class="msui-label">Name</span>
                  <TextInput size="sm" aria-label="Name" value={name} onInput={(v) => name.set(v)} />
                  <span class="msui-label">Quality</span>
                  <Select
                    size="sm"
                    over
                    aria-label="Quality"
                    options={[
                      ['auto', 'Auto'],
                      ['low', 'Low'],
                      ['medium', 'Medium'],
                      ['high', 'High'],
                    ]}
                    value={quality}
                    onChange={(v) => quality.set(v)}
                  />
                </div>
              </Popover>
            )}
          </Show>
        </div>
      );
    }),
};
