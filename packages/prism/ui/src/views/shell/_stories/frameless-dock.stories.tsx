/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { domView } from '../../../dom/view.js';
import type { ShellTab } from '../../../shell/shell.js';
import { mockBlankMode, ShellStory } from '../../../testing/shell.js';
import { story } from '../../../testing/story-root.js';

const meta: Meta = { title: 'Shell/Frameless dock', parameters: { layout: 'fullscreen' } };
export default meta;

const Text = 'MKLAVLTSEQ NCEDTRLWYH GGAPFQRSTV'.repeat(2);

/** A tab that draws on the canvas itself: text in the canvas layer's ink, no card around it. */
const tab = (frame: ShellTab['frame']): ShellTab => ({
  id: 'strip',
  title: 'Strip',
  icon: 'sequence',
  dock: 'top',
  frame,
  render: domView(() => (
    <div style="padding: var(--msp-space-3); font: 14px/1.6 var(--msp-font-mono, monospace); letter-spacing: 0.1em; color: var(--msp-color-text)">
      {Text}
    </div>
  )),
});

const Frameless = (background: string, frame: ShellTab['frame']): StoryObj => ({
  render: () =>
    story(() => (
      <div style="height: 100vh">
        <ShellStory
          loaded={true}
          background={background}
          mode={mockBlankMode()}
          addOn={({ shell }) => {
            queueMicrotask(() => shell.show('strip'));
            return { tabs: [tab(frame)] };
          }}
        />
      </div>
    )),
});

/** `frame: 'none'` on a light canvas: the text sits on the canvas, no glass card, shadow or hairline. */
export const Light: StoryObj = Frameless('#f4f5f7', 'none');

/** The same on a dark canvas: the ink follows the canvas background. */
export const Dark: StoryObj = Frameless('#1b1d22', 'none');

/** For comparison, the default (`frame: 'card'`): the glass card. */
export const Card: StoryObj = Frameless('#f4f5f7', 'card');
