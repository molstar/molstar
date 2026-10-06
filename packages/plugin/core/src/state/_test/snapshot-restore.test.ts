/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import type { PluginState } from '@molstar/plugin/state';
import { StructureRepresentation3D } from '@molstar/plugin/state/transforms/structure/representation';
import { StateTransform } from '@molstar/core/state';

const crambin = fs.readFileSync(path.resolve(__dirname, '../../../../../../examples/1crn.cif'), 'utf8');

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

/** The part of a snapshot that does not need a canvas. */
function takeSnapshot(plugin: PluginContext): PluginState.Snapshot {
  const snapshot = plugin.state.getSnapshot({
    camera: false,
    canvas3d: false,
    canvas3dContext: false,
    animation: false,
    startAnimation: false,
  } as any);
  // as saved to and loaded from a .molj file
  return JSON.parse(JSON.stringify(snapshot));
}

function cells(plugin: PluginContext) {
  return Array.from(plugin.state.data.cells.values()).filter((c) => c.transform.ref !== StateTransform.RootRef);
}

function representationTypes(plugin: PluginContext) {
  return cells(plugin)
    .filter((c) => c.transform.transformer === StructureRepresentation3D)
    .map((c) => (c.transform.params as any).type.name as string)
    .sort();
}

describe('restoring a snapshot into a plugin with the default spec', () => {
  it('rebuilds the data tree and the representations', async () => {
    const source = await createPlugin();
    await loadCrambin(source);
    const snapshot = takeSnapshot(source);
    const expectedTypes = representationTypes(source);
    expect(expectedTypes.length).toBeGreaterThan(0);
    expect(expectedTypes).toContain('cartoon');
    const expectedRefs = cells(source)
      .map((c) => c.transform.ref)
      .sort();
    source.dispose();

    const target = await createPlugin();
    expect(cells(target)).toHaveLength(0);
    await target.state.setSnapshot(snapshot);

    expect(
      cells(target)
        .map((c) => c.transform.ref)
        .sort(),
    ).toEqual(expectedRefs);
    expect(representationTypes(target)).toEqual(expectedTypes);
    expect(cells(target).every((c) => c.status === 'ok')).toBe(true);
    expect(target.managers.structure.hierarchy.current.structures).toHaveLength(1);
    target.dispose();
  });

  it('restores through the snapshot manager entry format as well', async () => {
    const source = await createPlugin();
    await loadCrambin(source);
    const snapshot = takeSnapshot(source);
    source.dispose();

    const target = await createPlugin();
    const state = {
      timestamp: 0,
      version: '6.0.0',
      entries: [{ snapshot, name: 'crambin' }],
      current: snapshot.id,
    };
    await target.managers.snapshot.setStateSnapshot(state as any);
    expect(target.managers.structure.hierarchy.current.structures).toHaveLength(1);
    target.dispose();
  });

  it('fails before changing anything when a transformer is not imported', async () => {
    const source = await createPlugin();
    await loadCrambin(source);
    const snapshot = takeSnapshot(source);
    source.dispose();

    const broken: PluginState.Snapshot = JSON.parse(JSON.stringify(snapshot));
    const transforms = broken.data!.tree.transforms as any[];
    const last = transforms[transforms.length - 1];
    (Array.isArray(last) ? last[0] : last).transformer = 'ms-plugin.not-imported-in-this-plugin';

    // the target already holds a state, which must survive the failed restore
    const target = await createPlugin();
    await loadCrambin(target);
    const refsBefore = cells(target).map((c) => c.transform.ref);
    const structuresBefore = target.managers.structure.hierarchy.current.structures.length;
    expect(refsBefore.length).toBeGreaterThan(0);

    await expect(target.state.setSnapshot(broken)).rejects.toThrow(
      /Snapshot uses transformers that are not available in this plugin: ms-plugin\.not-imported-in-this-plugin/,
    );
    expect(cells(target).map((c) => c.transform.ref)).toEqual(refsBefore);
    expect(target.managers.structure.hierarchy.current.structures).toHaveLength(structuresBefore);

    // the snapshot manager validates every entry before clearing itself
    const state = { timestamp: 0, version: '6.0.0', entries: [{ snapshot: broken }], current: broken.id };
    await expect(target.managers.snapshot.setStateSnapshot(state as any)).rejects.toThrow(/not available/);
    expect(target.managers.snapshot.state.entries.size).toBe(0);
    expect(cells(target).map((c) => c.transform.ref)).toEqual(refsBefore);
    target.dispose();
  });

  it('reports a representation that is not registered and renders the registry default', async () => {
    const source = await createPlugin();
    await loadCrambin(source);
    const snapshot = takeSnapshot(source);
    source.dispose();

    // a plugin that composed away cartoon
    const target = await createPlugin();
    const registry = target.representation.structure.registry;
    registry.remove(registry.get('cartoon'));
    expect(registry.has('cartoon')).toBe(false);
    const warn = jest.spyOn(target.log, 'warn').mockImplementation(() => {});
    await target.state.setSnapshot(snapshot);

    const messages = warn.mock.calls.map((c) => String(c[0]));
    expect(messages).toContain(
      "Structure representation 'cartoon' is not registered in this plugin; the registry default is used",
    );
    const types = representationTypes(target);
    expect(types).not.toContain('cartoon');
    expect(types).toContain(registry.default!.name);
    expect(cells(target).every((c) => c.status === 'ok')).toBe(true);
    target.dispose();
  });
});
