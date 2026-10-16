/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it, vi } from 'vitest';
import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginState } from '@molstar/plugin/state';
import { StateTransform } from '@molstar/core/state';

const crambin = fs.readFileSync(path.resolve(import.meta.dirname, '../../../../../../data/examples/1crn.cif'), 'utf8');

async function createPlugin() {
  const plugin = new PluginContext(DefaultPluginSpec());
  await plugin.init();
  return plugin;
}

async function loadCrambin(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: crambin, label: '1crn.cif' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
  await plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default');
}

const NoCanvas = {
  camera: false,
  canvas3d: false,
  canvas3dContext: false,
  animation: false,
  startAnimation: false,
} as any;

function cellCount(plugin: PluginContext) {
  return Array.from(plugin.state.data.cells.values()).filter((c) => c.transform.ref !== StateTransform.RootRef).length;
}

describe('PluginState.registerCustomState', () => {
  it('writes registered keys into snapshot.custom, whatever the snapshot params', async () => {
    const plugin = await createPlugin();
    expect(plugin.state.getSnapshot(NoCanvas).custom).toBeUndefined();

    plugin.state.registerCustomState('owner.a', { get: () => ({ version: 1, x: [1, 2] }), set: () => {} });
    plugin.state.registerCustomState('owner.b', { get: () => 'b', set: () => {} });
    const snapshot = plugin.state.getSnapshot({ ...NoCanvas, data: false });
    expect(snapshot.custom).toEqual({ 'owner.a': { version: 1, x: [1, 2] }, 'owner.b': 'b' });
    plugin.dispose();
  });

  it('calls set after the data tree is restored, with undefined for a missing key', async () => {
    const source = await createPlugin();
    await loadCrambin(source);
    const sourceCells = cellCount(source);
    expect(sourceCells).toBeGreaterThan(0);
    source.state.registerCustomState('owner.a', { get: () => ({ n: 1 }), set: () => {} });
    const snapshot: PluginState.Snapshot = JSON.parse(JSON.stringify(source.state.getSnapshot(NoCanvas)));
    source.dispose();

    const target = await createPlugin();
    const seen: { a: unknown[]; b: unknown[] } = { a: [], b: [] };
    let cellsWhenCalled = -1;
    target.state.registerCustomState('owner.a', {
      get: () => undefined,
      set: async (v) => {
        await Promise.resolve();
        cellsWhenCalled = cellCount(target);
        seen.a.push(v);
      },
    });
    target.state.registerCustomState('owner.b', { get: () => undefined, set: (v) => void seen.b.push(v) });
    await target.state.setSnapshot(snapshot);

    expect(seen.a).toEqual([{ n: 1 }]);
    expect(seen.b).toEqual([undefined]);
    expect(cellsWhenCalled).toBe(sourceCells);
    target.dispose();
  });

  it('ignores keys nobody registered and keeps them on the snapshot', async () => {
    const plugin = await createPlugin();
    const set = vi.fn();
    plugin.state.registerCustomState('owner.a', { get: () => 1, set });
    const snapshot = plugin.state.getSnapshot(NoCanvas);
    snapshot.custom = { ...snapshot.custom, 'someone.else': { keep: true } };
    await expect(plugin.state.setSnapshot(snapshot)).resolves.toBeUndefined();
    expect(set).toHaveBeenCalledWith(1);
    expect(snapshot.custom!['someone.else']).toEqual({ keep: true });
    plugin.dispose();
  });

  it('throws when a key is registered twice, and unregistering stops it', async () => {
    const plugin = await createPlugin();
    const set = vi.fn();
    const unregister = plugin.state.registerCustomState('owner.a', { get: () => 1, set });
    expect(() => plugin.state.registerCustomState('owner.a', { get: () => 2, set })).toThrow(/owner\.a/);

    unregister();
    expect(plugin.state.getSnapshot(NoCanvas).custom).toBeUndefined();
    await plugin.state.setSnapshot({ id: 'x' as any });
    expect(set).not.toHaveBeenCalled();

    // the key is free again
    plugin.state.registerCustomState('owner.a', { get: () => 3, set });
    expect(plugin.state.getSnapshot(NoCanvas).custom).toEqual({ 'owner.a': 3 });
    plugin.dispose();
  });

  it('keeps a failing set from breaking the restore', async () => {
    const plugin = await createPlugin();
    const error = vi.spyOn(plugin.log, 'error').mockImplementation(() => {});
    const set = vi.fn();
    plugin.state.registerCustomState('owner.a', {
      get: () => 1,
      set: () => {
        throw new Error('boom');
      },
    });
    plugin.state.registerCustomState('owner.b', { get: () => 2, set });
    await plugin.state.setSnapshot(plugin.state.getSnapshot(NoCanvas));
    expect(String(error.mock.calls[0][0])).toContain("'owner.a'");
    expect(set).toHaveBeenCalledWith(2);
    plugin.dispose();
  });

  it('does not apply custom state on transition frames', async () => {
    const plugin = await createPlugin();
    const set = vi.fn();
    plugin.state.registerCustomState('owner.a', { get: () => 1, set });
    const snapshot = plugin.state.getSnapshot(NoCanvas);
    const data = plugin.state.data.getSnapshot();
    snapshot.data = data;
    snapshot.transition = {
      frames: [
        { durationInMs: 100, data },
        { durationInMs: 100, data },
      ],
    };
    await plugin.state.setAnimationSnapshot(snapshot, 1);
    expect(set).not.toHaveBeenCalled();
    plugin.dispose();
  });

  it('survives the snapshot manager state and a JSON round trip', async () => {
    const source = await createPlugin();
    await loadCrambin(source);
    source.state.registerCustomState('owner.a', { get: () => ({ version: 1, edits: [['a', 2]] }), set: () => {} });
    const stateSnapshot = await source.managers.snapshot.getStateSnapshot({ params: NoCanvas });
    const json = JSON.parse(JSON.stringify(stateSnapshot));
    expect(json.entries[0].snapshot.custom).toEqual({ 'owner.a': { version: 1, edits: [['a', 2]] } });
    source.dispose();

    const target = await createPlugin();
    const set = vi.fn();
    target.state.registerCustomState('owner.a', { get: () => undefined, set });
    await target.managers.snapshot.setStateSnapshot(json);
    expect(set).toHaveBeenCalledWith({ version: 1, edits: [['a', 2]] });
    expect(target.managers.structure.hierarchy.current.structures).toHaveLength(1);
    target.dispose();
  });
});
