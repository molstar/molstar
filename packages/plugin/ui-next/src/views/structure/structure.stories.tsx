/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockSceneWithPicked } from '../../_stories/mocks/scene.js';
import { story } from '../../_stories/story-root.js';
import { ShellModel } from '../../shell/shell.js';
import { StructureTab } from './structure-tab.js';

const meta: Meta = { title: 'Views/Structure tab' };
export default meta;

const frame = 'width: 360px; height: 100vh; overflow: auto; background: var(--msp-color-panel); margin-left: auto;';

/** A representation picked: the inspector shows its parameters (essentials, then Advanced). */
export const Representation: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <StructureTab shell={new ShellModel()} scene={mockSceneWithPicked('polymer-cartoon')} />
      </div>
    )),
};

/** The structure picked. */
export const Structure: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <StructureTab shell={new ShellModel()} scene={mockSceneWithPicked('1tqn')} />
      </div>
    )),
};

/** Nothing picked: the tree only. */
export const TreeOnly: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <StructureTab shell={new ShellModel()} scene={mockSceneWithPicked('none')} />
      </div>
    )),
};
