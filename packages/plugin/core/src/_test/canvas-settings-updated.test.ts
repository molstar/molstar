/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { PluginCommands } from '@molstar/plugin/commands';
import { PluginContext } from '@molstar/plugin/context';

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

  it('a transition frame with canvas3d props emits', async () => {
    const { plugin, canvas3d, updated } = await createPlugin();
    await plugin.state.setAnimationSnapshot(
      { id: 'test', transition: { frames: [{ canvas3d: { props: {} } }] } } as any,
      0,
    );
    expect(canvas3d.setProps).toHaveBeenCalledTimes(1);
    expect(updated).toHaveBeenCalledTimes(1);
  });
});
