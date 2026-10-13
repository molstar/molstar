/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Asset } from '@molstar/core/util/assets';
import { Coordinates, Time } from '@molstar/model/model/structure';
import { PluginConfig } from '@molstar/plugin/config';
import { PluginContext } from '@molstar/plugin/context';
import { PluginStateObject as SO, PluginStateTransform } from '@molstar/plugin/state/objects';
import { DataFormatProvider } from '@molstar/plugin/state/formats/provider';
import { Mmcif } from '../../formats/trajectory/mmcif.js';
import { Pdb, PdbProvider } from '../../formats/trajectory/pdb.js';
import { Xyz, XyzProvider } from '../../formats/trajectory/xyz.js';
import { Cube, CubeProvider } from '../../formats/volume/cube.js';
import { Ccp4 } from '../../formats/volume/ccp4.js';
import { TrajectoryHierarchyPresetProvider } from '../../builder/structure/hierarchy-presets/types.js';
import { StructureRepresentationPresetProvider } from '../../builder/structure/representation-presets/types.js';
import { AddTrajectory, LoadTrajectory } from '../structure.js';

const HierarchyId = 'test-config-hierarchy';
const RepresentationId = 'test-config-representation';
const Missing = 'not-registered';
const MissingMessage = `Preset '${Missing}' is not registered in this plugin`;

const pdb = 'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  \n';
const cube = [
  'comment 1',
  'comment 2',
  '    1    0.000000    0.000000    0.000000',
  '    2    1.000000    0.000000    0.000000',
  '    2    0.000000    1.000000    0.000000',
  '    2    0.000000    0.000000    1.000000',
  '    1    0.000000    0.000000    0.000000    0.000000',
  '  0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8',
  '',
].join('\n');

const applied = { hierarchy: [] as unknown[], representation: [] as unknown[] };
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
  apply: async (_s, params) => {
    applied.representation.push(params);
    return {};
  },
});

/** A coordinates object of one frame of one atom. */
const TestCoordinates = PluginStateTransform.BuiltIn({
  name: 'test-config-presets-coordinates',
  display: 'Test coordinates',
  from: SO.Data.String,
  to: SO.Molecule.Coordinates,
})({
  apply() {
    const frame = {
      elementCount: 1,
      time: Time(0, 'ps'),
      x: [0],
      y: [0],
      z: [0],
      xyzOrdering: { isIdentity: true },
    };
    return new SO.Molecule.Coordinates(Coordinates.create([frame], Time(1, 'ps'), Time(0, 'ps')));
  },
});
const TestCoordinatesProvider = DataFormatProvider({
  name: 'test-config-presets-coordinates',
  label: 'Test coordinates',
  description: 'Test coordinates',
  category: 'Coordinates',
  stringExtensions: ['tcoords'],
  parse: async (plugin, data) => plugin.state.data.build().to(data).apply(TestCoordinates).commit(),
});

async function createPlugin() {
  const plugin = new PluginContext({ behaviors: [], registry: [] });
  await plugin.init();
  plugin.dataFormats.clear();
  plugin.register([
    Pdb,
    Xyz,
    Cube,
    Ccp4,
    Mmcif,
    {
      formats: [TestCoordinatesProvider],
      structure: { presets: { hierarchy: [hierarchy], representation: [representation] } },
    },
  ]);
  jest
    .spyOn(plugin.builders.data, 'download')
    .mockImplementation(({ label }) =>
      plugin.builders.data.rawData({ data: label === 'cube' ? cube : pdb, label }, { state: { isGhost: true } }),
    );
  applied.hierarchy.length = 0;
  applied.representation.length = 0;
  return plugin;
}

/** Runs `f` with the transaction error log captured. */
async function captureErrors(f: () => Promise<unknown>) {
  const error = jest.spyOn(console, 'error').mockImplementation(() => {});
  try {
    await f();
    return error.mock.calls.map((c) => String(c[0]));
  } finally {
    error.mockRestore();
  }
}

describe('trajectory format visuals', () => {
  it('apply the configured hierarchy preset', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, HierarchyId);
    const data = await plugin.builders.data.download({ url: 'a.pdb' } as any);
    const parsed = await PdbProvider.parse(plugin, data);
    await PdbProvider.visuals!(plugin, parsed);
    expect(applied.hierarchy.length).toBe(1);

    await XyzProvider.visuals!(plugin, parsed);
    expect(applied.hierarchy.length).toBe(2);
  });

  it('fail when the configured hierarchy preset is not registered', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultHierarchyPreset, Missing);
    const data = await plugin.builders.data.download({ url: 'a.pdb' } as any);
    const parsed = await PdbProvider.parse(plugin, data);
    expect(() => PdbProvider.visuals!(plugin, parsed)).toThrow(MissingMessage);
  });

  it('use the default hierarchy preset when nothing is configured', async () => {
    const plugin = await createPlugin();
    const data = await plugin.builders.data.download({ url: 'a.pdb' } as any);
    const parsed = await PdbProvider.parse(plugin, data);
    const spy = jest.spyOn(plugin.builders.structure.hierarchy, 'applyPreset').mockResolvedValue({} as any);
    await PdbProvider.visuals!(plugin, parsed);
    expect(spy).toHaveBeenCalledWith(expect.anything(), 'preset-trajectory-default');
  });
});

describe('Cube format visuals', () => {
  async function load(plugin: PluginContext) {
    const data = await plugin.builders.data.download({ label: 'cube' } as any);
    const parsed = await CubeProvider.parse(plugin, data);
    return CubeProvider.visuals!(plugin, parsed);
  }

  it('applies the configured representation preset to the structure', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, RepresentationId);
    const visuals = await load(plugin);
    expect(applied.representation.length).toBe(1);
    // one isosurface, shown by the entry's providers
    expect(visuals.length).toBe(1);
  });

  it('fails when the configured representation preset is not registered', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, Missing);
    await expect(load(plugin)).rejects.toThrow(MissingMessage);
  });
});

describe('structure hierarchy manager', () => {
  async function createStructure(plugin: PluginContext) {
    const data = await plugin.builders.data.download({ url: 'a.pdb' } as any);
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    return plugin.builders.structure.createStructure(model);
  }

  it('applies the configured representation preset when the structure type changes', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, RepresentationId);
    const structure = await createStructure(plugin);
    await plugin.managers.structure.hierarchy.updateStructure({ cell: structure.cell! } as any, {
      name: 'model',
      params: {},
    });
    expect(applied.representation.length).toBe(1);
  });

  it('fails when the configured representation preset is not registered', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, Missing);
    const structure = await createStructure(plugin);
    const errors = await captureErrors(() =>
      plugin.managers.structure.hierarchy.updateStructure({ cell: structure.cell! } as any, {
        name: 'model',
        params: {},
      }),
    );
    expect(errors.join('\n')).toContain(MissingMessage);
  });
});

describe('trajectory actions', () => {
  async function createModelAndCoordinates(plugin: PluginContext) {
    const data = await plugin.builders.data.download({ url: 'a.pdb' } as any);
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
    const model = await plugin.builders.structure.createModel(trajectory);
    const text = await plugin.builders.data.rawData({ data: 'x', label: 'coords' }, { state: { isGhost: true } });
    const coordinates = await plugin.state.data.build().to(text).apply(TestCoordinates).commit();
    return { model, coordinates };
  }

  function addTrajectory(plugin: PluginContext, model: string, coordinates: string) {
    return plugin.runTask(plugin.state.data.applyAction(AddTrajectory, { model, coordinates }));
  }

  it('AddTrajectory applies the configured representation preset', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, RepresentationId);
    const { model, coordinates } = await createModelAndCoordinates(plugin);
    await addTrajectory(plugin, model.ref, coordinates.ref);
    expect(applied.representation.length).toBe(1);
  });

  it('AddTrajectory fails when the configured representation preset is not registered', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, Missing);
    const { model, coordinates } = await createModelAndCoordinates(plugin);
    const errors = await captureErrors(() => addTrajectory(plugin, model.ref, coordinates.ref));
    expect(errors.join('\n')).toContain(MissingMessage);
  });

  function loadTrajectory(plugin: PluginContext) {
    const params = {
      source: {
        name: 'url' as const,
        params: {
          model: { url: Asset.Url('https://example.org/a.pdb'), format: 'pdb', isBinary: false },
          coordinates: { url: Asset.Url('https://example.org/a.tcoords'), format: TestCoordinatesProvider.name },
        },
      },
    };
    return plugin.runTask(plugin.state.data.applyAction(LoadTrajectory, params as any));
  }

  it('LoadTrajectory applies the configured representation preset', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, RepresentationId);
    await loadTrajectory(plugin);
    expect(applied.representation.length).toBe(1);
  });

  it('LoadTrajectory fails when the configured representation preset is not registered', async () => {
    const plugin = await createPlugin();
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, Missing);
    const logged: string[] = [];
    plugin.events.log.subscribe((e) => {
      if (e.type === 'error') logged.push(e.message);
    });
    const errors = await captureErrors(() => loadTrajectory(plugin));
    // the action logs the failure
    expect(errors.join('\n')).toContain(MissingMessage);
    expect(logged).toEqual(['Error loading trajectory']);
  });
});
