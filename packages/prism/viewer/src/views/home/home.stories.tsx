/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockHome } from '../../testing/home.js';
import { story } from '@molstar/prism-ui/testing/index';
import { HomeTab } from './home-tab.js';

const meta: Meta = { title: 'Views/Home tab' };
export default meta;

const frame = 'width: 360px; height: 100vh; overflow: auto; background: var(--msp-color-panel); margin-left: auto;';

/** The root's actions as forms (the first open, Open Files among them), and the command menu hint. */
export const Default: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <HomeTab home={mockHome((text) => console.info(text))} />
      </div>
    )),
};

/** With remote states (a feature flag in the spec, off by default): sessions on a state server. */
export const WithRemoteStates: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <HomeTab home={mockHome((text) => console.info(text), { remoteStates: true })} />
      </div>
    )),
};
