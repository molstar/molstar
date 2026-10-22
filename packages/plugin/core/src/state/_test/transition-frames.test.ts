/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { PluginState } from '@molstar/plugin/state';

/** A transition whose frames last 100, 1000 and 100 ms: frame 1 is a hold */
function snapshot(loop = false): PluginState.Snapshot {
  const data = { tree: { transforms: [] } } as any;
  const frames = [100, 1000, 100].map((durationInMs) => ({ durationInMs, data }));
  return { id: 'a' as any, transition: { frames, loop } };
}

describe('the frame of a transition at a time', () => {
  it('is the frame whose slot holds the time', () => {
    const s = snapshot();
    const at = (t: number) => PluginState.getStateTransitionFrameIndex(s, t);
    expect([0, 50, 99.9].map(at)).toEqual([0, 0, 0]);
    expect([100, 150, 600, 1099.9].map(at)).toEqual([1, 1, 1, 1]);
    expect([1100, 1150, 1199.9].map(at)).toEqual([2, 2, 2]);
  });

  it('is the last frame after the end, and the first before the start', () => {
    const s = snapshot();
    expect(PluginState.getStateTransitionFrameIndex(s, 1200)).toBe(2);
    expect(PluginState.getStateTransitionFrameIndex(s, 5000)).toBe(2);
    expect(PluginState.getStateTransitionFrameIndex(s, -10)).toBe(0);
  });

  it('wraps around when the transition loops', () => {
    const s = snapshot(true);
    expect([1200, 1250, 1300, 2350].map((t) => PluginState.getStateTransitionFrameIndex(s, t))).toEqual([0, 0, 1, 2]);
  });

  it('starts where the frames before it end', () => {
    const s = snapshot();
    expect([0, 1, 2].map((i) => PluginState.getStateTransitionFrameTime(s, i))).toEqual([0, 100, 1100]);
  });

  it('is undefined without frames', () => {
    expect(PluginState.getStateTransitionFrameIndex({ id: 'b' as any }, 0)).toBeUndefined();
  });
});
