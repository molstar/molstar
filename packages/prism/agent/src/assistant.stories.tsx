/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { story } from '@molstar/prism-ui/testing/index';
import { ViewerShellStory } from '@molstar/prism-viewer/testing/index';
import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockAgentAddOn } from './testing/index.js';

const meta: Meta = { title: 'Shell/Assistant', parameters: { layout: 'fullscreen' } };
export default meta;

/**
 * The agent module included and set up (a mock agent stands in for any-harness): ask from the command menu (⌘K, type
 * a request, Ask the assistant; “color it red” opens the filled sentence, “show the active site”, “load 1w0e and
 * superpose it”, “try something broken” run turns), watch the turn card, Allow or Deny a download, Stop, Undo turn;
 * History at the rail's foot; the Assistant tab with the model and keys.
 */
export const SetUp: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} addOn={mockAgentAddOn({ configured: true })} />
      </div>
    )),
};

/** The agent module included but no model yet: the Assistant tab asks for a key; the Ask row opens it. */
export const NotSetUp: StoryObj = {
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ViewerShellStory loaded={true} addOn={mockAgentAddOn()} />
      </div>
    )),
};
