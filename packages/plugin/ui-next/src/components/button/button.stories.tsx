/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { map } from 'rxjs';
import { story } from '../../_stories/story-root.js';
import { observableState } from '../../reactive/observable-state.js';
import { Button } from './button.js';

const meta: Meta = { title: 'Components/Button' };
export default meta;

/** Placeholder glyph until the icon set lands. */
const Dot = () => (
  <svg viewBox="0 0 24 24" fill="currentColor">
    <circle cx="12" cy="12" r="6" />
  </svg>
);

const row = 'display: flex; gap: 8px; align-items: center; flex-wrap: wrap;';

export const Variants: StoryObj = {
  render: () =>
    story(() => {
      const pressed = observableState(false);
      return (
        <div style="display: grid; gap: 12px; padding: 16px">
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
            <Button icon={<Dot />}>With icon</Button>
            <Button icon={<Dot />} title="Icon only" />
            <Button icon={<Dot />} size="sm" title="Small icon only" variant="ghost" />
          </div>
          <div style={row}>
            <Button pressed={pressed} onClick={() => pressed.set(!pressed.value)}>
              Toggle
            </Button>
            <span>pressed: {pressed.asReadonly.pipe(map(String))}</span>
          </div>
        </div>
      );
    }),
};
