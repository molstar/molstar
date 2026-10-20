/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import { PluginContext } from '@molstar/plugin/context';
import { ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { PluginSpec } from '@molstar/plugin/spec';
import { SnapshotControls } from '../../../behavior/dynamic/state.js';
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
    const navigator = { applyKey: vi.fn(), applyNext: vi.fn(), applyId: vi.fn() };
    const remove = manager.setNavigator(navigator);
    const set = vi.spyOn(plugin.state, 'setSnapshot');
    manager.applyKey('a');
    manager.applyNext(-1);
    const first = manager.state.entries.get(0)!.snapshot.id;
    manager.applyId(first);
    expect(navigator.applyKey).toHaveBeenCalledWith('a');
    expect(navigator.applyNext).toHaveBeenCalledWith(-1);
    expect(navigator.applyId).toHaveBeenCalledWith(first);
    expect(set).not.toHaveBeenCalled();

    remove();
    manager.applyKey('a');
    expect(set).toHaveBeenCalledTimes(1);
    expect(manager.getEntry(manager.state.current!)?.key).toBe('a');
    await manager.applyId(manager.state.entries.get(1)!.snapshot.id);
    expect(set).toHaveBeenCalledTimes(2);
    expect(manager.getEntry(manager.state.current!)?.key).toBe('b');
    expect(manager.applyId('no-such-id')).toBeUndefined();
    plugin.dispose();
  });

  it("the keyboard's first and last go through applyId; with no snapshot they do nothing", async () => {
    const plugin = new PluginContext({ behaviors: [PluginSpec.Behavior(SnapshotControls)] });
    await plugin.init();
    (plugin as unknown as { canvas3d: object }).canvas3d = {};
    const manager = plugin.managers.snapshot;
    const navigator = { applyKey: vi.fn(), applyNext: vi.fn(), applyId: vi.fn() };
    manager.setNavigator(navigator);
    const press = (code: string) =>
      plugin.behaviors.interaction.keyReleased.next({
        code,
        key: code,
        modifiers: ModifiersKeys.create({ control: true }),
      } as never);
    expect(() => press('ArrowUp')).not.toThrow();
    expect(navigator.applyId).not.toHaveBeenCalled();
    for (const key of ['a', 'b', 'c']) {
      manager.add(PluginStateSnapshotManager.Entry(plugin.state.getSnapshot({ camera: false }), { key, name: key }));
    }
    const ids = manager.state.entries.map((e) => e.snapshot.id).toArray();
    press('ArrowUp');
    press('ArrowDown');
    expect(navigator.applyId.mock.calls).toEqual([[ids[0]], [ids[2]]]);
    // The stand-in canvas has nothing to dispose.
    (plugin as unknown as { canvas3d: undefined }).canvas3d = undefined;
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
