/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PrismAppStory, story } from '@molstar/prism-ui/testing/index';
import { mockViewerMode } from '@molstar/prism-viewer/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockStoryMode } from '../testing/index.js';

const meta: Meta = { title: 'Shell/Story mode', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * An app with two modes, started in the story mode (a first sketch: the Story tab with the scenes and the caption, the
 * scenes playing on the canvas, story commands in ⌘K). "Switch mode…" in the command menu shows the choice over the
 * app, saying what goes; choosing the viewer starts a new app in it, Cancel goes back.
 */
export const StartInStory: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <PrismAppStory loaded={true} modes={[mockStoryMode(), mockViewerMode()]} start="story" />
      </div>
    )),
};

/** The viewer mode first, the story mode a switch away. */
export const StartInViewer: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <PrismAppStory loaded={true} modes={[mockViewerMode(), mockStoryMode()]} start="viewer" />
      </div>
    )),
};

/**
 * The first open: the splash asks which mode before any app exists; choosing one creates the app in it. With nothing
 * loaded, Home offers the modes again (no loss: it starts the other at once); later, "Switch mode…" asks first.
 */
export const FirstOpen: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <PrismAppStory loaded={false} modes={[mockViewerMode(), mockStoryMode()]} />
      </div>
    )),
};
