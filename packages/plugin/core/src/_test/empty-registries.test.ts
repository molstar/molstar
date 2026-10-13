/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';

const Scopes = ['structure', 'volume', 'particles'] as const;

const names = (list: readonly { name: string }[]) => list.map((e) => e.name);

function actionIds(manager: unknown) {
  const { actions, fromTypeIndex } = manager as { actions: Map<string, unknown>; fromTypeIndex: Map<unknown, unknown> };
  return { count: actions.size, types: fromTypeIndex.size };
}

/** What a plugin lists in every registry that `PluginSpec.registry` can fill. */
function contents(plugin: PluginContext) {
  const r = plugin.representation;
  return {
    representations: Object.fromEntries(Scopes.map((s) => [s, names(r[s].registry.list)])),
    themes: Object.fromEntries(
      Scopes.map((s) => [
        s,
        {
          color: names(r[s].themes.colorThemeRegistry.list),
          size: names(r[s].themes.sizeThemeRegistry.list),
        },
      ]),
    ),
    formats: names(plugin.dataFormats.list),
    hierarchyPresets: plugin.builders.structure.hierarchy.providers.map((p) => p.id),
    representationPresets: plugin.builders.structure.representation.providers.map((p) => p.id),
    selectionQueries: plugin.query.structure.registry.list.length,
    markdownExtensions: names((plugin.managers.markdownExtensions as any).extension),
    dragAndDrop: plugin.managers.dragAndDrop.list().map((e) => e.name),
    animations: names(plugin.managers.animation.animations),
    dataActions: actionIds(plugin.state.data.actions).count,
    lociLabelProviders: plugin.managers.lociLabels.providers.length,
  };
}

describe('a plugin with an empty spec', () => {
  it('starts with empty registries and managers', async () => {
    const plugin = new PluginContext({ behaviors: [] });
    await plugin.init();

    for (const scope of Scopes) {
      const { registry, themes } = plugin.representation[scope];
      expect(registry.list).toEqual([]);
      expect(registry.default).toBeUndefined();
      expect(themes.colorThemeRegistry.list).toEqual([]);
      expect(themes.sizeThemeRegistry.list).toEqual([]);
    }
    expect(plugin.dataFormats.list).toEqual([]);
    expect(plugin.dataFormats.extensions.size).toBe(0);
    expect(plugin.builders.structure.hierarchy.providers).toEqual([]);
    expect(plugin.builders.structure.representation.providers).toEqual([]);
    expect(plugin.query.structure.registry.list).toEqual([]);
    expect(plugin.query.structure.registry.options).toEqual([]);
    expect((plugin.managers.markdownExtensions as any).extension).toEqual([]);
    expect(plugin.managers.dragAndDrop.list()).toEqual([]);
    expect(plugin.managers.animation.animations).toEqual([]);
    expect(plugin.managers.animation.isEmpty).toBe(true);
    expect(plugin.managers.lociLabels.providers).toEqual([]);
    expect(actionIds(plugin.state.data.actions)).toEqual({ count: 0, types: 0 });
    expect(actionIds(plugin.state.behaviors.actions)).toEqual({ count: 0, types: 0 });
    plugin.dispose();
  });

  it('does not use a built-in provider when a name is requested', async () => {
    const plugin = new PluginContext({ behaviors: [] });
    await plugin.init();
    expect(plugin.dataFormats.has('mmcif')).toBe(false);
    expect(() => plugin.dataFormats.get('mmcif')).toThrow(/not registered/);
    expect(plugin.builders.structure.hierarchy.has('default')).toBe(false);
    expect(plugin.builders.structure.representation.has('auto')).toBe(false);
    expect(() => plugin.builders.structure.hierarchy.applyPreset({} as any, 'default')).toThrow(/not registered/);
    plugin.dispose();
  });
});

describe('the default spec', () => {
  it('registers the 5.x provider sets in the 5.x order', async () => {
    const plugin = new PluginContext(DefaultPluginSpec());
    await plugin.init();
    const c = contents(plugin);

    expect(c.representations).toEqual({
      structure: [
        'cartoon',
        'backbone',
        'ball-and-stick',
        'blob-surface',
        'carbohydrate',
        'ellipsoid',
        'gaussian-surface',
        'gaussian-volume',
        'label',
        'line',
        'molecular-surface',
        'orientation',
        'plane',
        'point',
        'putty',
        'spacefill',
        'polyhedron',
        'interactions',
        'cross-link-restraint',
      ],
      volume: ['direct-volume', 'dot', 'isosurface', 'segment', 'slice', 'streamlines'],
      particles: ['spacefill', 'orientation', 'fibers', 'target'],
    });

    // 5.x built-in themes (37 color in the volume and particle scopes, 41 in the structure scope with the themes the
    // behaviors add), plus the external structure and volume themes in every scope
    expect(c.themes.structure.color).toHaveLength(43);
    expect(c.themes.volume.color).toHaveLength(39);
    expect(c.themes.particles.color).toHaveLength(39);
    for (const scope of Scopes) {
      expect(c.themes[scope].color).toEqual(
        expect.arrayContaining(['uniform', 'external-structure', 'external-volume']),
      );
      expect(c.themes[scope].size).toHaveLength(6);
    }

    expect(c.formats).toEqual([
      'ccp4',
      'dsn6',
      'cube',
      'dx',
      'dscif',
      'segcif',
      'sfcif',
      'mtz',
      'psf',
      'prmtop',
      'top',
      'dcd',
      'xtc',
      'trr',
      'nctraj',
      'lammpstrj',
      'ply',
      'obj',
      'vtp',
      'relion_star_particles',
      'dynamo_tbl_particles',
      'cryoet_ndjson_particles',
      'artiatomi_em_particles',
      'mmcif_particles',
      'simularium_particles',
      'mmcif',
      'cifCore',
      'pdb',
      'pdbqt',
      'pqr',
      'gro',
      'xyz',
      'lammps_data',
      'lammps_traj_data',
      'mol',
      'sdf',
      'mol2',
    ]);
    expect(c.hierarchyPresets).toEqual([
      'preset-trajectory-default',
      'preset-trajectory-all-models',
      'preset-trajectory-unitcell',
      'preset-trajectory-supercell',
      'preset-trajectory-crystal-contacts',
    ]);
    expect(c.representationPresets).toEqual([
      'preset-structure-representation-empty',
      'preset-structure-representation-auto',
      'preset-structure-representation-atomic-detail',
      'preset-structure-representation-polymer-cartoon',
      'preset-structure-representation-polymer-and-ligand',
      'preset-structure-representation-protein-and-nucleic',
      'preset-structure-representation-coarse-surface',
      'preset-structure-representation-illustrative',
      'preset-structure-representation-molecular-surface',
      'preset-structure-representation-auto-lod',
      'preset-structure-representation-mesoscale',
    ]);
    expect(c.selectionQueries).toBe(67);
    expect(c.markdownExtensions).toEqual([
      'center-camera',
      'apply-snapshot',
      'next-snapshot',
      'focus-refs',
      'highlight-refs',
      'query',
      'play-audio',
      'toggle-audio',
      'pause-audio',
      'stop-audio',
      'dispose-audio',
      'play-transition',
      'play-snapshots',
      'stop-animation',
    ]);
    // the open-anything handler, which 5.x hard-coded in the manager, is an entry now
    expect(c.dragAndDrop).toEqual(['open-files']);
    expect(c.animations).toEqual([
      'built-in.animate-model-index',
      'built-in.animate-particle-trajectory',
      'built-in.animate-camera-spin',
      'built-in.animate-camera-rock',
      'built-in.animate-state-snapshots',
      'built-in.animate-state-snapshot-transition',
      'built-in.animate-assembly-unwind',
      'built-in.animate-structure-spin',
      'built-in.animate-state-interpolation',
      'built-in.animate-time',
    ]);
    expect(c.dataActions).toBe(58);
    expect(c.lociLabelProviders).toBe(5);
    plugin.dispose();
  });
});
