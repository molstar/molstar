/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { PluginStateSnapshotManager } from '../snapshots.js';

async function createPlugin() {
  const plugin = new PluginContext({ behaviors: [] });
  await plugin.init();
  const manager = plugin.managers.snapshot;
  for (const key of ['a', 'b']) {
    manager.add(PluginStateSnapshotManager.Entry(plugin.state.getSnapshot({ camera: false }), { key, name: key }));
  }
  return { plugin, manager };
}

describe("the snapshot manager's navigator", () => {
  it('takes applyKey and applyNext while set; the manager applies them again once it is removed', async () => {
    const { plugin, manager } = await createPlugin();
    const navigator = { applyKey: vi.fn(), applyNext: vi.fn() };
    const remove = manager.setNavigator(navigator);
    const set = vi.spyOn(plugin.state, 'setSnapshot');
    manager.applyKey('a');
    manager.applyNext(-1);
    expect(navigator.applyKey).toHaveBeenCalledWith('a');
    expect(navigator.applyNext).toHaveBeenCalledWith(-1);
    expect(set).not.toHaveBeenCalled();

    remove();
    manager.applyKey('a');
    expect(set).toHaveBeenCalledTimes(1);
    expect(manager.getEntry(manager.state.current!)?.key).toBe('a');
    plugin.dispose();
  });

  it('clears the current snapshot, saying so once', async () => {
    const { plugin, manager } = await createPlugin();
    expect(manager.state.current).toBeDefined();
    const changed = vi.fn();
    manager.events.changed.subscribe(changed);
    manager.clearCurrent();
    manager.clearCurrent();
    expect(manager.state.current).toBeUndefined();
    expect(changed).toHaveBeenCalledTimes(1);
    plugin.dispose();
  });
});
