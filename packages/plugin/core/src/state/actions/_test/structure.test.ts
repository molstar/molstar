/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginConfig } from '@molstar/plugin/config';
import { PluginContext } from '@molstar/plugin/context';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { DefaultFormats } from '@molstar/plugin/default-registry';
import { Ccp4 } from '../../formats/volume/ccp4.js';
import { Mmcif } from '../../formats/trajectory/mmcif.js';
import { Mol } from '../../formats/trajectory/mol.js';
import { Pdb } from '../../formats/trajectory/pdb.js';
import { Sdf } from '../../formats/trajectory/sdf.js';
import { Xyz } from '../../formats/trajectory/xyz.js';
import { TrajectoryHierarchyPresetProvider } from '../../builder/structure/hierarchy-presets/types.js';
import { StructureRepresentationPresetProvider } from '../../builder/structure/representation-presets/types.js';
import { EmptyPresetEntry } from '../../builder/structure/representation-presets/empty.js';
import { DownloadStructure } from '../structure.js';

const AllSources = ['pdb', 'pdb-ihm', 'swissmodel', 'alphafolddb', 'modelarchive', 'pubchem', 'url'];

/** A plugin with the given format entries only. */
async function createPlugin(...formats: PluginRegistryEntry[]) {
  const plugin = new PluginContext({ behaviors: [], registry: formats });
  await plugin.init();
  return plugin;
}

function getParams(plugin: PluginContext) {
  const params = DownloadStructure.definition.params!(void 0 as any, plugin) as { source: PD.Mapped<any> };
  const { source } = params;
  const sources = source.select.options.map((o) => o[0] as string);
  const map = (name: string) => source.map(name) as PD.Group<any>;
  return {
    source,
    sources,
    map,
    urlFormat: sources.includes('url') ? (map('url').params.format as PD.Select<string>) : undefined,
    asTrajectory: sources[0] ? (map(sources[0]).params.options?.params.asTrajectory as PD.Base<boolean>) : undefined,
  };
}

describe('DownloadStructure params', () => {
  it('offers every source and the trajectory formats of the registry by default', async () => {
    const plugin = await createPlugin(DefaultFormats);
    const { sources, urlFormat, source, asTrajectory } = getParams(plugin);
    expect(sources).toEqual(AllSources);
    expect(source.defaultValue.name).toBe('pdb');
    expect(urlFormat!.defaultValue).toBe('mmcif');
    expect(urlFormat!.options.map((o) => o[0])).toEqual(
      plugin.dataFormats.list.filter((f) => f.provider.category === 'Trajectory').map((f) => f.name),
    );
    expect(urlFormat!.options.map((o) => o[0])).toEqual([
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
    expect(asTrajectory!.isHidden).toBeFalsy();
  });

  it('offers only the sources whose formats are registered', async () => {
    const mmcif = getParams(await createPlugin(Mmcif));
    expect(mmcif.sources).toEqual(['pdb', 'pdb-ihm', 'alphafolddb', 'modelarchive', 'url']);
    expect(mmcif.urlFormat!.options.map((o) => o[0])).toEqual(['mmcif']);
    expect(mmcif.urlFormat!.defaultValue).toBe('mmcif');
    expect(mmcif.asTrajectory!.isHidden).toBeFalsy();

    const pdb = getParams(await createPlugin(Pdb));
    expect(pdb.sources).toEqual(['swissmodel', 'url']);
    expect(pdb.source.defaultValue.name).toBe('swissmodel');
    expect(pdb.urlFormat!.options.map((o) => o[0])).toEqual(['pdb']);
    expect(pdb.urlFormat!.defaultValue).toBe('pdb');

    const mol = getParams(await createPlugin(Mol, Sdf));
    expect(mol.sources).toEqual(['pubchem', 'url']);
    expect(mol.urlFormat!.options.map((o) => o[0])).toEqual(['mol', 'sdf']);
    expect(mol.urlFormat!.defaultValue).toBe('mol');

    const other = getParams(await createPlugin(Xyz));
    expect(other.sources).toEqual(['url']);
    expect(other.urlFormat!.defaultValue).toBe('xyz');
  });

  it('hides the multi-source blob option unless mmCIF is registered', async () => {
    const plugin = await createPlugin(Pdb);
    const { map } = getParams(plugin);
    const options = map('swissmodel').params.options.params.asTrajectory as PD.Base<boolean>;
    expect(options.isHidden).toBe(true);
  });

  it('offers nothing when no trajectory format is registered', async () => {
    const { sources, source } = getParams(await createPlugin(Ccp4));
    expect(sources).toEqual(['']);
    expect(source.defaultValue.name).toBe('');
  });

  it('defaults the representation preset to the configured one', async () => {
    const plugin = await createPlugin(Mmcif, EmptyPresetEntry);
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, 'preset-structure-representation-empty');
    const { map } = getParams(plugin);
    expect(map('pdb').params.options.params.representation.defaultValue).toBe('preset-structure-representation-empty');
  });
});

describe('DownloadStructure presets', () => {
  const HierarchyId = 'test-config-hierarchy';
  const RepresentationId = 'test-config-representation';
  const pdb = ['ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  '].join('\n');

  const applied: { hierarchy: any[] } = { hierarchy: [] };
  const hierarchy = TrajectoryHierarchyPresetProvider({
    id: HierarchyId,
    display: { name: 'Test hierarchy' },
    apply: async (_t, params) => {
      applied.hierarchy.push(params);
      return {};
    },
  });
  const representation = StructureRepresentationPresetProvider({
    id: RepresentationId,
    display: { name: 'Test representation' },
    apply: async () => ({}),
  });

  async function createDownloadPlugin() {
    const plugin = await createPlugin(Pdb, EmptyPresetEntry);
    plugin.register({ structure: { presets: { hierarchy: [hierarchy], representation: [representation] } } });
    jest
      .spyOn(plugin.builders.data, 'download')
      .mockImplementation(() =>
        plugin.builders.data.rawData({ data: pdb, label: 'a.pdb' }, { state: { isGhost: true } }),
      );
    applied.hierarchy.length = 0;
    return plugin;
  }

  function run(plugin: PluginContext) {
    const params = DownloadStructure.createDefaultParams(void 0 as any, plugin) as any;
    params.source = {
      name: 'url',
      params: { ...params.source.params, url: Asset.Url('https://example.org/a.pdb'), format: 'pdb' },
    };
    return plugin.runTask(plugin.state.data.applyAction(DownloadStructure, params));
  }

  it('applies the configured hierarchy preset with the configured representation preset', async () => {
    const plugin = await createDownloadPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, HierarchyId);
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, RepresentationId);
    await run(plugin);
    expect(applied.hierarchy.length).toBe(1);
    expect(applied.hierarchy[0]).toMatchObject({ representationPreset: RepresentationId, showUnitcell: true });
  });

  it('does not show the unit cell for the empty representation preset', async () => {
    const plugin = await createDownloadPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, HierarchyId);
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, 'preset-structure-representation-empty');
    await run(plugin);
    expect(applied.hierarchy[0]).toMatchObject({
      representationPreset: 'preset-structure-representation-empty',
      showUnitcell: false,
    });
  });

  it('fails when the configured hierarchy preset is not registered', async () => {
    const plugin = await createDownloadPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, 'not-registered');
    // the action's transaction reverts and logs the error
    const error = jest.spyOn(console, 'error').mockImplementation(() => {});
    await run(plugin);
    expect(error.mock.calls.map((c) => String(c[0]))).toEqual([
      "Error: Preset 'not-registered' is not registered in this plugin",
    ]);
    error.mockRestore();
    expect(applied.hierarchy.length).toBe(0);
    // the downloaded data was reverted
    expect(plugin.state.data.cells.size).toBe(1);
  });

  it('uses the default hierarchy preset of the registry by default', async () => {
    const plugin = await createDownloadPlugin();
    const spy = jest.spyOn(plugin.builders.structure.hierarchy, 'applyPreset').mockResolvedValue({} as any);
    await run(plugin);
    expect(spy).toHaveBeenCalledWith(expect.anything(), 'preset-trajectory-default', expect.anything());
  });
});
