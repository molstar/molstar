/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockHome } from '../../_stories/mocks/home.js';
import { story } from '../../_stories/story-root.js';
import { HomeTab } from './home-tab.js';

const meta: Meta = { title: 'Views/Home tab' };
export default meta;

const frame = 'width: 360px; height: 100vh; overflow: auto; background: var(--msp-color-panel); margin-left: auto;';

/** The root's actions as forms (the first open), recent entries, open files, and the command menu hint. */
export const Default: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <HomeTab home={mockHome((text) => console.info(text))} />
      </div>
    )),
};
