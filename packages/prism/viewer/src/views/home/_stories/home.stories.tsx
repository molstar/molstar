/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { story } from '@molstar/prism-ui/testing/index';
import { RealViewerStory } from '../../../testing/real-viewer.js';

// Home in the real viewer: what loads data (the loaders' actions as forms, Open files), opened at start while nothing is
// loaded (`ui.homeWhenEmpty`).

const meta: Meta = { title: 'Viewer/Home', parameters: { layout: 'fullscreen' } };
export default meta;

/** The root's actions as forms (the first open, Open files among them), and the command menu hint. */
export const Default: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory spec={(base) => ({ ...base, options: { ...base.options, ui: { homeWhenEmpty: true } } })} />
    )),
};

/** With remote states (on in the spec, off by default): sessions on a state server. Network. */
export const WithRemoteStates: StoryObj = {
  render: () =>
    story(() => (
      <RealViewerStory
        spec={(base) => ({
          ...base,
          options: { ...base.options, ui: { homeWhenEmpty: true, showRemoteStates: true } },
        })}
      />
    )),
};
