/** @jsxImportSource @molstar/plugin-ui-next/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockSnapshots, mockTrajectory } from '../../_stories/mocks/playback.js';
import { story } from '../../_stories/story-root.js';
import { PlaybackCard } from './playback.js';

const meta: Meta = { title: 'Canvas/Playback' };
export default meta;

const frame = 'width: 320px; padding: 24px; background: var(--msp-color-bg);';
const post = (text: string) => console.info(text);

/** Snapshots (dots for a few) and the trajectory below a hairline. */
export const SnapshotsAndTrajectory: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <PlaybackCard snapshots={mockSnapshots(post)} trajectory={mockTrajectory()} />
      </div>
    )),
};

/** Many snapshots: the position is "n / N", opening a filterable list of titles. */
export const ManySnapshots: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <PlaybackCard snapshots={mockSnapshots(post, 32)} />
      </div>
    )),
};

/** Only a trajectory: the card is one row. */
export const TrajectoryOnly: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <PlaybackCard trajectory={mockTrajectory(50)} />
      </div>
    )),
};
