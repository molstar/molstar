/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '../../../testing/story-root.js';
import { observableState } from '../../../reactive/observable-state.js';
import { Icon } from '../../icon/icon.js';
import { Button } from '../button.js';

import { derive } from '../../../reactive/state.js';
const meta: Meta = { title: 'Components/Button' };
export default meta;

const row = 'display: flex; gap: var(--msp-space-3); align-items: center; flex-wrap: wrap;';

export const Variants: StoryObj = {
  render: () =>
    story(() => {
      const pressed = observableState(false);
      return (
        <div style="display: grid; gap: var(--msp-space-4); padding: var(--msp-space-5)">
          <div style={row}>
            <Button>Default</Button>
            <Button variant="primary">Apply</Button>
            <Button variant="ghost">Ghost</Button>
            <Button variant="danger">Remove</Button>
            <Button disabled>Disabled</Button>
          </div>
          <div style={row}>
            <Button size="sm">Small</Button>
            <Button size="sm" variant="primary">
              Small primary
            </Button>
            <Button size="xs">Extra small</Button>
            <Button icon={<Icon name="refresh-cw" />}>With icon</Button>
            <Button icon={<Icon name="plus" />} title="Icon only" />
            <Button icon={<Icon name="trash-2" />} size="sm" title="Small icon only" variant="ghost" />
            <Button icon={<Icon name="refresh-cw" />} size="xs" title="Extra small icon only" variant="ghost" />
          </div>
          <div style="width: 300px">
            <Button block variant="primary">
              Block
            </Button>
          </div>
          <div style={row}>
            <Button pressed={pressed} onClick={() => pressed.set(!pressed.value)}>
              Toggle
            </Button>
            <span>pressed: {derive(pressed.asReadonly, String)}</span>
          </div>
        </div>
      );
    }),
};
