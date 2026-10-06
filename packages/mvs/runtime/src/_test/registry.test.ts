/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import { setFSModule } from '@molstar/core/util/data-source';
import { PluginContext } from '@molstar/plugin/context';
import { PluginSpec, type PluginRegistryEntry } from '@molstar/plugin/spec';
import { MVSData } from '@molstar/mvs-builder/mvs-data';
import {
  MVSRepresentationParams,
  MVSVolumeRepresentationParams,
} from '@molstar/mvs-builder/tree/mvs/mvs-tree-representations';
import { MolViewSpec, createMVSRegistryEntry } from '@molstar/mvs/behavior';
import { MVS_BEHAVIOR_ID } from '@molstar/mvs/behavior-id';
import { loadMVS } from '@molstar/mvs/load';
import { MVSRuntimeRegistry } from '@molstar/mvs/registry';

setFSModule(fs);

const LocalStructure = `file://${process.cwd()}/examples/1cbs_full.bcif`;

type Scope = 'structure' | 'volume';
type Names = { representations: Set<string>; color: Set<string>; size: Set<string> };
const emptyNames = (): Names => ({ representations: new Set(), color: new Set(), size: new Set() });

const TransformerScopes: Record<string, Scope> = {
  'ms-plugin.structure-representation-3d': 'structure',
  'ms-plugin.volume-representation-3d': 'volume',
};

const nameOf = (v: any): string | undefined => (v && typeof v.name === 'string' ? v.name : undefined);
const names = (providers: readonly { name: string }[] | undefined) => (providers ?? []).map((p) => p.name);

/** Collects the representation and theme names the MVS loader wrote into the data trees of all snapshots. */
function collectEmittedNames(plugin: PluginContext) {
  const emitted: Record<Scope, Names> = { structure: emptyNames(), volume: emptyNames() };
  const entries = plugin.managers.snapshot.state.entries.toArray();
  for (const entry of entries) {
    for (const t of entry.snapshot.data?.tree.transforms ?? []) {
      const scope = TransformerScopes[t.transformer];
      if (!scope) continue;
      const params = t.params as any;
      const type = nameOf(params?.type);
      const color = nameOf(params?.colorTheme);
      const size = nameOf(params?.sizeTheme);
      if (type) emitted[scope].representations.add(type);
      if (color) emitted[scope].color.add(color);
      if (size) emitted[scope].size.add(size);
    }
  }
  return { emitted, snapshots: entries.length };
}

/** An MVS document with one snapshot per way MVS names a representation or a theme. */
function createCoverageDocument() {
  const snapshots: ReturnType<typeof MVSData.stateToStates>['snapshots'] = [];
  const add = (build: (root: ReturnType<typeof MVSData.createBuilder>) => void) => {
    const builder = MVSData.createBuilder();
    build(builder);
    snapshots.push(...MVSData.stateToStates(builder.getState()).snapshots);
  };
  const annotation = { uri: 'annotations.json', format: 'json', schema: 'all_atomic' } as const;

  // The first snapshot is the one `loadMVS` applies, so it loads from the local file.
  add((b) => {
    const struct = b.download({ url: LocalStructure }).parse({ format: 'bcif' }).modelStructure();
    struct.component({ selector: 'polymer' }).representation({ type: 'cartoon' }).color({ color: 'green' });
  });

  // Every representation type of the schema and its variants, with uniform and split colors.
  add((b) => {
    const struct = b.download({ url: 'structure.bcif' }).parse({ format: 'bcif' }).modelStructure();
    const comp = struct.component({ selector: 'all' });
    for (const type of Object.keys(MVSRepresentationParams.cases) as (keyof typeof MVSRepresentationParams.cases)[]) {
      comp.representation({ type }).color({ color: 'red' });
    }
    comp.representation({ type: 'surface', surface_type: 'gaussian' }).color({ color: 'red/blue' });
    comp.representation({ type: 'putty', size_theme: 'uncertainty' }).color({ color: 'red' });
    comp.representation({ type: 'putty', size_theme: 'uniform' }).color({ color: 'red' });
  });

  // Several color layers, annotation colors, labels, and the non-covalent interactions extension.
  add((b) => {
    const struct = b.download({ url: 'structure.bcif' }).parse({ format: 'bcif' }).modelStructure();
    struct
      .component({ selector: 'polymer' })
      .representation({ type: 'cartoon' })
      .color({ color: 'green' })
      .color({ color: 'blue', selector: { label_asym_id: 'A' } });
    struct.component({ selector: 'polymer' }).representation({ type: 'spacefill' }).colorFromUri(annotation);
    struct
      .component({ selector: 'polymer' })
      .representation({ type: 'line' })
      .colorFromSource({ schema: 'all_atomic' });
    struct.labelFromUri(annotation);
    struct.labelFromSource({ schema: 'all_atomic' });
    struct.component({ selector: 'ligand' }).label({ text: 'Ligand' });
    struct.component({ selector: 'ligand', custom: { molstar_show_non_covalent_interactions: true } });
  });

  // Volume representations.
  add((b) => {
    const volume = b.download({ url: 'map.ccp4' }).parse({ format: 'map' }).volume();
    for (const type of Object.keys(MVSVolumeRepresentationParams.cases)) {
      if (type === 'isosurface') volume.representation({ type }).color({ color: 'blue' });
      else if (type === 'grid_slice') volume.representation({ type, dimension: 'x' }).color({ color: 'blue' });
      else throw new Error(`The coverage document does not cover the volume representation type '${type}'`);
    }
  });

  return MVSData.createMultistate(snapshots);
}

async function createPlugin(registry: PluginRegistryEntry[]) {
  const plugin = new PluginContext({ registry, behaviors: [PluginSpec.Behavior(MolViewSpec)] });
  await plugin.init();
  return plugin;
}

describe('MVSRuntimeRegistry', () => {
  it('lists the representations of the MVS nodes and interactions, in the structure and volume scopes only', () => {
    expect(Object.keys(MVSRuntimeRegistry).sort()).toEqual(['structure', 'volume']);
    expect(names(MVSRuntimeRegistry.structure?.representations).sort()).toEqual(
      [
        'backbone',
        'ball-and-stick',
        'carbohydrate',
        'cartoon',
        'gaussian-surface',
        'interactions',
        'line',
        'molecular-surface',
        'putty',
        'spacefill',
      ].sort(),
    );
    expect(names(MVSRuntimeRegistry.volume?.representations).sort()).toEqual(['isosurface', 'slice']);
  });

  it('registers every representation and theme name that MVS can emit', async () => {
    // The plugin has the behavior, which registers what MVS itself defines (its labels and color themes).
    const plugin = await createPlugin([MVSRuntimeRegistry]);
    const doc = createCoverageDocument();
    await loadMVS(plugin, doc);
    expect(plugin.errorContext.get('mvs')).toEqual([]);
    const { emitted, snapshots } = collectEmittedNames(plugin);
    expect(snapshots).toBe(doc.snapshots.length);

    // What the loader emits; a new name appears here and in the registry (or in the behavior entry).
    expect([...emitted.structure.representations].sort()).toEqual(
      [
        'backbone',
        'ball-and-stick',
        'carbohydrate',
        'cartoon',
        'mvs-custom-label',
        'gaussian-surface',
        'interactions',
        'line',
        'molecular-surface',
        'mvs-annotation-label',
        'putty',
        'spacefill',
      ].sort(),
    );
    expect([...emitted.structure.color].sort()).toEqual(
      ['element-symbol', 'interaction-type', 'mvs-annotation', 'mvs-multilayer', 'mvs-split-uniform', 'uniform'].sort(),
    );
    expect([...emitted.structure.size].sort()).toEqual(['physical', 'uncertainty', 'uniform']);
    expect([...emitted.volume.representations].sort()).toEqual(['isosurface', 'slice']);
    expect([...emitted.volume.color]).toEqual(['uniform']);

    const { structure, volume } = plugin.representation;
    for (const name of emitted.structure.representations) expect(structure.registry.has(name)).toBe(true);
    for (const name of emitted.structure.color) expect(structure.themes.colorThemeRegistry.has(name)).toBe(true);
    for (const name of emitted.structure.size) expect(structure.themes.sizeThemeRegistry.has(name)).toBe(true);
    for (const name of emitted.volume.representations) expect(volume.registry.has(name)).toBe(true);
    for (const name of emitted.volume.color) expect(volume.themes.colorThemeRegistry.has(name)).toBe(true);
  });

  it('lists exactly what MVS names: each provider is emitted by MVS or is the default theme of a listed representation', async () => {
    const plugin = await createPlugin([MVSRuntimeRegistry]);
    await loadMVS(plugin, createCoverageDocument());
    const { emitted } = collectEmittedNames(plugin);
    // The providers the behavior registers itself.
    const own = createMVSRegistryEntry(plugin.representation.structure.themes.colorThemeRegistry).structure!;
    const ownRepresentations = new Set(names(own.representations));
    const ownColorThemes = new Set(names(own.themes?.color));

    for (const scope of ['structure', 'volume'] as const) {
      const entry = MVSRuntimeRegistry[scope]!;
      expect(names(entry.representations).sort()).toEqual(
        [...emitted[scope].representations].filter((n) => !(scope === 'structure' && ownRepresentations.has(n))).sort(),
      );

      const defaultColors = entry.representations!.map((r) => r.defaultColorTheme.name);
      const defaultSizes = entry.representations!.map((r) => r.defaultSizeTheme.name);
      const colors = new Set([
        ...[...emitted[scope].color].filter((n) => !(scope === 'structure' && ownColorThemes.has(n))),
        ...defaultColors,
      ]);
      const sizes = new Set([...emitted[scope].size, ...defaultSizes]);
      expect(names(entry.themes?.color).sort()).toEqual([...colors].sort());
      expect(names(entry.themes?.size).sort()).toEqual([...sizes].sort());
    }
  });

  it('brings the default themes of every representation it lists', () => {
    for (const scope of ['structure', 'volume'] as const) {
      const entry = MVSRuntimeRegistry[scope]!;
      const color = new Set(names(entry.themes?.color));
      const size = new Set(names(entry.themes?.size));
      for (const r of entry.representations!) {
        expect(color.has(r.defaultColorTheme.name)).toBe(true);
        expect(size.has(r.defaultSizeTheme.name)).toBe(true);
      }
    }
  });
});

describe('MolViewSpec behavior entry', () => {
  const state = (plugin: PluginContext) => ({
    representations: plugin.representation.structure.registry.list.map((r) => r.name),
    colorThemes: plugin.representation.structure.themes.colorThemeRegistry.list.map((r) => r.name).sort(),
    formats: plugin.dataFormats.list.map((f) => f.name),
    dragAndDrop: plugin.managers.dragAndDrop.list().map((e) => e.name),
    lociLabels: plugin.managers.lociLabels.providers.length,
  });

  it('holds the provider part of the behavior', async () => {
    const plugin = await createPlugin([]);
    const entry = createMVSRegistryEntry(plugin.representation.structure.themes.colorThemeRegistry);
    expect(Object.keys(entry).sort()).toEqual(['actions', 'dragAndDrop', 'formats', 'lociLabels', 'structure']);
    expect(names(entry.structure?.representations)).toEqual(['mvs-custom-label', 'mvs-annotation-label']);
    expect(names(entry.structure?.themes?.color)).toEqual(['mvs-split-uniform', 'mvs-annotation', 'mvs-multilayer']);
    expect(names(entry.formats)).toEqual(['MVSJ', 'MVSX']);
    expect(entry.lociLabels?.length).toBe(2);
    expect(entry.dragAndDrop?.map((d) => d.name)).toEqual(['mvs-mvsj-mvsx']);
    expect(entry.actions?.length).toBe(1);
  });

  it('registers the entry when the behavior is on and removes exactly it when the behavior is turned off', async () => {
    const plugin = await createPlugin([]);
    expect(state(plugin)).toEqual({
      representations: ['mvs-custom-label', 'mvs-annotation-label'],
      colorThemes: ['mvs-annotation', 'mvs-multilayer', 'mvs-split-uniform'],
      formats: ['MVSJ', 'MVSX'],
      dragAndDrop: ['mvs-mvsj-mvsx'],
      lociLabels: 2,
    });

    // A provider the plugin also registers for itself survives the behavior being turned off.
    const unregisterShared = plugin.register({
      structure: {
        representations: createMVSRegistryEntry(plugin.representation.structure.themes.colorThemeRegistry).structure!
          .representations,
      },
    });

    await plugin.state.behaviors.build().delete(MVS_BEHAVIOR_ID).commit();
    expect(plugin.state.hasBehavior(MVS_BEHAVIOR_ID)).toBe(false);
    expect(state(plugin)).toEqual({
      representations: ['mvs-custom-label', 'mvs-annotation-label'],
      colorThemes: [],
      formats: [],
      dragAndDrop: [],
      lociLabels: 0,
    });
    unregisterShared();
    expect(state(plugin).representations).toEqual([]);
  });
});
