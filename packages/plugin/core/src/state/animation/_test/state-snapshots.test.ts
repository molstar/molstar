/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Camera } from '@molstar/graphics/canvas3d/camera';
import {
  AnimateStateSnapshots,
  AnimateStateSnapshotTransition,
} from '@molstar/plugin/state/animation/built-in/state-snapshots';
import { addTransition, createPlaybackPlugin, play, rawDataTree } from './harness.js';

/** Three looks of a 4-node tree, as frames of 100 ms each: a, a, b, b, c, c (holds share the frame's data) */
async function transitionPlugin(durationInMs?: number) {
  const env = await createPlaybackPlugin();
  const { plugin } = env;
  const a = await rawDataTree(plugin, 4, 'a');
  const b = await rawDataTree(plugin, 4, 'b');
  const c = await rawDataTree(plugin, 4, 'c');
  const frames = [a, a, b, b, c, c].map((data) => ({ durationInMs: 100, data }));
  const camera = new Camera().getSnapshot();
  addTransition(plugin, frames, {
    camera: { current: camera, transitionStyle: 'animate', transitionDurationInMs: 250 },
    durationInMs,
  });
  return env;
}

describe('a playing snapshot transition', () => {
  it('applies each frame once, and moves the camera once', async () => {
    const { plugin, canvas3d, counts } = await transitionPlugin();
    const updated = counts.updated;
    await plugin.managers.animation.play(AnimateStateSnapshotTransition, {});
    const setup = counts.setSnapshot;
    const ticks = await play(plugin, 700);
    expect(ticks).toBe(37);
    // a frame index change applies its frame: frame 0 once with the camera (setup applied it without), then 2, 3, 4, 5
    // (frame 1 is skipped: B2)
    expect(counts.setSnapshot - setup).toBe(5);
    // the camera of the snapshot is set with frame 0, once: a reset per tick would restart the transition
    expect(canvas3d.requestCameraReset).toHaveBeenCalledTimes(1);
    expect(canvas3d.requestCameraReset.mock.calls[0][0]).toMatchObject({ durationMs: 250 });
    // setup goes back from look c to a, then b and c: the frames that change the tree update its four nodes
    expect(counts.updated - updated).toBe(12);
    expect(counts.recreated).toBe(0);
  });

  it('stays on the last frame when it ends', async () => {
    const { plugin, counts } = await transitionPlugin();
    await plugin.managers.animation.play(AnimateStateSnapshotTransition, {});
    await play(plugin, 2000);
    expect(plugin.managers.animation.isAnimating).toBe(false);
    expect(
      [...plugin.state.data.cells.values()].filter((c) => c.transform.ref !== '-=root=-').map((c) => c.obj?.label),
    ).toEqual(['c0', 'c1', 'c2', 'c3']);
    expect(counts.recreated).toBe(0);
  });
});

describe('playing the snapshots', () => {
  it('applies each frame of a transition once', async () => {
    const { plugin, canvas3d, counts } = await transitionPlugin(600);
    await plugin.managers.animation.play(AnimateStateSnapshots, {});
    const setup = counts.setSnapshot;
    await play(plugin, 590);
    // frames 2..5 (frame 1 is skipped: B2); frame 0 of the first snapshot is the snapshot itself, which setup applied
    expect(counts.setSnapshot - setup).toBe(4);
    expect(canvas3d.requestCameraReset).toHaveBeenCalledTimes(1);
  });
});
