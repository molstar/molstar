/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockSnapshots, mockTrajectory } from '../../../testing/playback.js';
import { story } from '../../../testing/story-root.js';
import { SnapshotCard, TrajectoryBar } from '../playback.js';

const meta: Meta = { title: 'Canvas/Playback' };
export default meta;

const frame = 'width: 320px; padding: var(--msp-space-6); background: var(--msp-color-bg);';
const post = (text: string) => console.info(text);

/** Two different things, two surfaces: the session's snapshots, and one structure's models. */
export const SnapshotsAndTrajectory: StoryObj = {
  render: () =>
    story(() => (
      <div style={`${frame} display: grid; gap: var(--msp-space-3);`}>
        <SnapshotCard model={mockSnapshots(post)} />
        <TrajectoryBar model={mockTrajectory()} />
      </div>
    )),
};

/** Many snapshots: the position is "n / N", opening a filterable list of titles. */
export const ManySnapshots: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <SnapshotCard model={mockSnapshots(post, 32)} />
      </div>
    )),
};

/** A trajectory: the bar names its structure. */
export const Trajectory: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <TrajectoryBar model={mockTrajectory(50, '6M0J')} />
      </div>
    )),
};
