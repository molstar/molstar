/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ShellStory, story } from '@molstar/prism-ui/testing/index';
import { mockViewerMode } from '@molstar/prism-viewer/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockStoryMode } from './testing/index.js';

const meta: Meta = { title: 'Shell/Story mode', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * The app with two modes: story mode (a first sketch: the Story tab with the scenes and the caption, the scenes playing
 * on the canvas, story commands in ⌘K) and the viewer mode (Mol*'s own UI). The rail's switch (above the theme) or the command menu
 * (“Switch to …”) swaps them: the frame stays (status, Log, Settings, the theme), everything a mode brought goes.
 */
export const SwitchModes: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={true} modes={[mockStoryMode(), mockViewerMode()]} />
      </div>
    )),
};

/** The viewer mode first, the story mode a switch away. */
export const StartInViewer: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={true} modes={[mockViewerMode(), mockStoryMode()]} />
      </div>
    )),
};

/**
 * The first open: the app asks which mode to start in, over the empty canvas. Choosing one starts it; afterwards the
 * rail's switch changes modes, and Home offers the choice again while nothing is loaded (Clear empties the session).
 */
export const FirstOpen: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory loaded={false} chooseMode={true} modes={[mockViewerMode(), mockStoryMode()]} />
      </div>
    )),
};
