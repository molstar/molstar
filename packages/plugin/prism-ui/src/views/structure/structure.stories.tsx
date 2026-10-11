/** @jsxImportSource @molstar/plugin-prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import {
  type MockScenePreset,
  mockComponentOptions,
  mockPresets,
  mockQuickStyles,
  mockSceneWithPicked,
} from '../../_stories/mocks/scene.js';
import { story } from '../../_stories/story-root.js';
import { ShellModel } from '../../shell/shell.js';
import { StructureTab } from './structure-tab.js';

const meta: Meta = { title: 'Views/Structure tab' };
export default meta;

const frame = 'width: 360px; height: 100vh; overflow: auto; background: var(--msp-color-panel); margin-left: auto;';

function Tab(props: { pick: string; preset?: MockScenePreset; stream?: boolean }) {
  const shell = new ShellModel();
  const scene = mockSceneWithPicked(props.pick, { preset: props.preset });
  const picked = scene.picked.value;
  if (props.stream && picked) scene.streaming(picked).set('enabled', true);
  return (
    <div style={frame}>
      <StructureTab
        shell={shell}
        scene={scene}
        presets={mockPresets()}
        quickStyles={mockQuickStyles()}
        options={mockComponentOptions()}
      />
    </div>
  );
}

/** A representation picked: the inspector shows its parameters (essentials, then Advanced). */
export const Representation: StoryObj = { render: () => story(() => <Tab pick="polymer-cartoon" />) };

/** The structure picked: facts, its parameters, then extension sections (density streaming; no symmetry for C1). */
export const Structure: StoryObj = { render: () => story(() => <Tab pick="1tqn" />) };

/** Density streaming enabled from the structure's section: the channels appear under the structure. */
export const StreamingDensity: StoryObj = { render: () => story(() => <Tab pick="1tqn" stream />) };

/** Several structures, collapsed but the first. */
export const ManyStructures: StoryObj = { render: () => story(() => <Tab pick="none" preset="many" />) };

/** An EMDB map next to its fitted model; the map picked: its facts and the isosurface below it. */
export const EmdbMap: StoryObj = { render: () => story(() => <Tab pick="emd-2984-iso" preset="em" />) };

/** The fitted model of a symmetric assembly: assembly symmetry shows. */
export const EmdbModel: StoryObj = { render: () => story(() => <Tab pick="5a1a" preset="em" />) };

/** Nothing picked: the tree only. */
export const TreeOnly: StoryObj = { render: () => story(() => <Tab pick="none" />) };
