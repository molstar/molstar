/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { AnimateStateInterpolation } from '@molstar/plugin/state/animation/built-in/state-interpolation';
import { PluginStateSnapshotManager } from '@molstar/plugin/state/manager/snapshots';
import { createPlaybackPlugin, play, rawDataTree } from './harness.js';

/** Two snapshots of a 4-node tree whose labels differ; raw data has no `interpolate`, so it steps at the middle */
async function twoSnapshots() {
  const env = await createPlaybackPlugin();
  const { plugin } = env;
  const b = await rawDataTree(plugin, 4, 'b');
  const a = await rawDataTree(plugin, 4, 'a');
  for (const [id, data] of [
    ['a', a],
    ['b', b],
  ] as const) {
    plugin.managers.snapshot.add(PluginStateSnapshotManager.Entry({ id: id as any, data }, { name: id }));
  }
  return env;
}

describe('the state interpolation animation', () => {
  it('updates a node when its interpolated params change, not on every tick', async () => {
    const { plugin, counts } = await twoSnapshots();
    const updated = counts.updated;
    const busy: boolean[] = [];
    plugin.behaviors.state.isUpdating.subscribe((b) => busy.push(b));
    busy.length = 0;
    await plugin.managers.animation.play(AnimateStateInterpolation, { transtionDurationInMs: 1000 });
    // 0 to 950 ms: the source (a, as shown) until the middle, then the target (b): the four nodes change once, in one
    // update; the other ticks commit nothing
    const ticks = await play(plugin, 950);
    expect(ticks).toBe(58);
    expect(counts.updated - updated).toBe(4);
    expect(busy.filter((b) => b)).toHaveLength(1);
    expect(plugin.state.data.cells.get('node0')!.obj!.label).toBe('b0');
  });
});
