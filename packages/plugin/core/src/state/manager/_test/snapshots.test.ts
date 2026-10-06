/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { StateTransformer } from '@molstar/core/state';
import type { PluginContext } from '@molstar/plugin/context';
import type { PluginState } from '@molstar/plugin/state';
import { PluginStateSnapshotManager } from '../snapshots.js';

const Valid = StateTransformer.ROOT.id;

function snapshot(id: string, transformer: string): PluginState.Snapshot {
  return {
    id: id as any,
    data: { tree: { transforms: [{ parent: 'ref', transformer, params: void 0, ref: 'ref', version: '1' }] } },
  };
}

function stateSnapshot(...entries: PluginState.Snapshot[]): PluginStateSnapshotManager.StateSnapshot {
  return {
    timestamp: 1,
    version: 'test',
    current: void 0,
    playback: { isPlaying: false, nextSnapshotDelayInMs: 1000 },
    entries: entries.map((s) => PluginStateSnapshotManager.Entry(s, {})),
  };
}

describe('PluginStateSnapshotManager.setStateSnapshot', () => {
  function create() {
    const setSnapshot = jest.fn(async (_snapshot: PluginState.Snapshot) => {});
    const plugin = { state: { setSnapshot }, managers: { asset: { delete: jest.fn() } } } as unknown as PluginContext;
    return { manager: new PluginStateSnapshotManager(plugin), setSnapshot };
  }

  it('loads a valid multi-entry state', async () => {
    const { manager, setSnapshot } = create();
    const a = snapshot('a', Valid);
    const b = snapshot('b', Valid);
    expect(await manager.setStateSnapshot(stateSnapshot(a, b))).toBe(a);
    expect(manager.state.entries.size).toBe(2);
    expect(setSnapshot).toHaveBeenCalledWith(a);
    manager.dispose();
  });

  it('leaves the manager unchanged when any entry names an unregistered transformer', async () => {
    const { manager, setSnapshot } = create();
    await manager.setStateSnapshot(stateSnapshot(snapshot('old-1', Valid), snapshot('old-2', Valid)));
    setSnapshot.mockClear();
    const before = manager.state;
    const changed = jest.fn();
    manager.events.changed.subscribe(changed);

    const next = stateSnapshot(snapshot('new-1', Valid), snapshot('new-2', 'test.missing'));
    await expect(manager.setStateSnapshot(next)).rejects.toThrow(
      'Snapshot uses transformers that are not available in this plugin: test.missing.',
    );

    expect(manager.state).toBe(before);
    expect(manager.state.entries.map((e) => e.snapshot.id).toArray()).toEqual(['old-1', 'old-2']);
    expect(manager.getEntry('old-1')).toBeDefined();
    expect(manager.getEntry('new-1')).toBeUndefined();
    expect(changed).not.toHaveBeenCalled();
    expect(setSnapshot).not.toHaveBeenCalled();
    manager.dispose();
  });
});
