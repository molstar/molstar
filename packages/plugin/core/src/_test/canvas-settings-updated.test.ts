/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginContext } from '@molstar/plugin/context';
import { AnimateStateSnapshotTransition } from '@molstar/plugin/state/animation/built-in/state-snapshots';

/** A plugin without WebGL: `canvas3d` and `canvas3dContext` are fakes that record `setProps` calls. */
async function createPlugin() {
  const plugin = new PluginContext({ behaviors: [] });
  await plugin.init();
  const canvas3d = { props: {}, setProps: vi.fn() };
  const canvas3dContext = { props: {}, setProps: vi.fn() };
  (plugin as any).canvas3d = canvas3d;
  (plugin as any).canvas3dContext = canvas3dContext;
  const updated = vi.fn();
  plugin.events.canvas3d.settingsUpdated.subscribe(updated);
  return { plugin, canvas3d, canvas3dContext, updated };
}

describe('canvas props emit settingsUpdated', () => {
  it('setCanvas3DProps sets the props and emits once', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    plugin.setCanvas3DProps({ checkeredTransparentBackground: true });
    expect(canvas3d.setProps).toHaveBeenCalledTimes(1);
    expect(updated).toHaveBeenCalledTimes(1);
  });

  it('setCanvas3DContextProps sets the props and emits once', async () => {
    const { plugin, canvas3dContext, updated } = await createPlugin();
    plugin.setCanvas3DContextProps({ transparency: 'blended' });
    expect(canvas3dContext.setProps).toHaveBeenCalledWith({ transparency: 'blended' });
    expect(updated).toHaveBeenCalledTimes(1);
  });

  it('neither emits without a canvas', async () => {
    const plugin = new PluginContext({ behaviors: [] });
    await plugin.init();
    const updated = vi.fn();
    plugin.events.canvas3d.settingsUpdated.subscribe(updated);
    plugin.setCanvas3DProps({});
    plugin.setCanvas3DContextProps({});
    expect(updated).not.toHaveBeenCalled();
  });

  it('clear(true) resets the viewport settings and emits once', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    await plugin.clear(true);
    expect(canvas3d.setProps).toHaveBeenCalledTimes(1);
    expect(updated).toHaveBeenCalledTimes(1);
  });

  it('clear() leaves the canvas alone', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    await plugin.clear();
    expect(canvas3d.setProps).not.toHaveBeenCalled();
    expect(updated).not.toHaveBeenCalled();
  });

  it('the SetSettings command emits once', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    await PluginCommands.Canvas3D.SetSettings(plugin, { settings: { checkeredTransparentBackground: true } });
    expect(canvas3d.setProps).toHaveBeenCalledTimes(1);
    expect(updated).toHaveBeenCalledTimes(1);
  });

  it('snapshot restore of canvas3d and canvas3dContext props emits', async () => {
    const { plugin, canvas3d, canvas3dContext, updated } = await createPlugin();
    await plugin.state.setSnapshot({
      id: 'test',
      canvas3d: { props: {} as any },
      canvas3dContext: { props: {} as any },
    } as any);
    expect(canvas3d.setProps).toHaveBeenCalledTimes(1);
    expect(canvas3dContext.setProps).toHaveBeenCalledTimes(1);
    // one emit for the SetSettings command, one for the context props
    expect(updated).toHaveBeenCalledTimes(2);
  });

  it('scrubbed transition frames with canvas3d props emit once they settle', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    vi.useFakeTimers();
    try {
      const snapshot = { id: 'test', transition: { frames: [{ canvas3d: { props: {} } }] } } as any;
      await plugin.state.setAnimationSnapshot(snapshot, 0);
      await plugin.state.setAnimationSnapshot(snapshot, 0);
      expect(canvas3d.setProps).toHaveBeenCalledTimes(2);
      expect(updated).not.toHaveBeenCalled();
      vi.runAllTimers();
      expect(updated).toHaveBeenCalledTimes(1);
    } finally {
      vi.useRealTimers();
    }
  });

  it('a playing transition sets the frames silently and emits once on teardown', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    const frame = { durationInMs: 10, canvas3d: { props: {} } };
    const snapshot = { id: 'test', canvas3d: { props: {} }, transition: { frames: [frame, frame, frame] } } as any;
    const state = { totalDuration: 30, snapshots: [{ snapshot }], currentIndex: 0, isInitial: true };
    const ctx = { params: {}, plugin };
    for (const current of [5, 15]) {
      await AnimateStateSnapshotTransition.apply(state, { lastApplied: 0, current }, ctx);
    }
    expect(canvas3d.setProps).toHaveBeenCalledTimes(2);
    expect(updated).not.toHaveBeenCalled();
    await AnimateStateSnapshotTransition.teardown!({}, state, plugin);
    expect(updated).toHaveBeenCalledTimes(1);
    // nothing new to say
    await AnimateStateSnapshotTransition.teardown!({}, state, plugin);
    expect(updated).toHaveBeenCalledTimes(1);
  });
});
