/** @jsxImportSource @molstar/prism-ui/dom */

/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Meta, StoryObj } from '@storybook/web-components-vite';
import { mockPeriodicTrajectory, mockSnapshots, mockTrajectory, mockTransition } from '../../../testing/playback.js';
import { story } from '../../../testing/story-root.js';
import { PlaybackDeck, SnapshotCard, SnapshotRow, TrajectoryRow, TransitionRow } from '../playback.js';

// The playback views over stand-in models, for their states one by one (the viewer's own: Viewer › Playback): the
// canvas's playback card and its rows, and the Story card's snapshot.

const meta: Meta = { title: 'Components/Playback' };
export default meta;

const frame = 'width: 320px; padding: var(--msp-space-6); background: var(--msp-color-bg);';
const post = (text: string) => console.info(text);

/** The playback card: a row each for the snapshots, the current one's transition and a trajectory's models. */
export const Card: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <PlaybackDeck>
          <SnapshotRow model={mockSnapshots(post)} />
          <TransitionRow model={mockTransition()} />
          <TrajectoryRow model={mockTrajectory(50, '6M0J')} />
        </PlaybackDeck>
      </div>
    )),
};

/** A simulation in a periodic box: its frames, and the box menu (molecules whole, in the box, no jumps, as read). */
export const PeriodicBox: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <PlaybackDeck>
          <TrajectoryRow model={mockPeriodicTrajectory()} />
        </PlaybackDeck>
      </div>
    )),
};

/** Many snapshots: the position is "n / N", opening a filterable list of titles. */
export const ManySnapshots: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <PlaybackDeck>
          <SnapshotRow model={mockSnapshots(post, 32)} />
        </PlaybackDeck>
      </div>
    )),
};

/** The Story card's snapshot: its title, its position as dots, its description (markdown, links to others). */
export const Snapshot: StoryObj = {
  render: () =>
    story(() => (
      <div style={frame}>
        <SnapshotCard model={mockSnapshots(post)} />
      </div>
    )),
};
