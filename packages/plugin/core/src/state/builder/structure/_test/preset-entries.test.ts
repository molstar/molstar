/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import * as fs from 'node:fs';
import * as path from 'node:path';
import { PluginContext } from '@molstar/plugin/context';
import { PluginConfig } from '@molstar/plugin/config';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import '@molstar/plugin/state/transforms/catalog';
import { Structure } from '@molstar/model/model/structure';
import { DefaultPresets } from '@molstar/plugin/default-registry';
import { PresetStructureRepresentations } from '../representation-presets/catalog.js';
import { PresetTrajectoryHierarchy } from '../hierarchy-presets/catalog.js';
import { AutoPresetEntry } from '../representation-presets/auto.js';
import { AtomicDetailPresetEntry } from '../representation-presets/atomic-detail.js';
import { AutoLodPresetEntry } from '../representation-presets/auto-lod.js';
import { CoarseSurfacePresetEntry } from '../representation-presets/coarse-surface.js';
import { EmptyPresetEntry } from '../representation-presets/empty.js';
import { IllustrativePresetEntry } from '../representation-presets/illustrative.js';
import { MesoscalePresetEntry } from '../representation-presets/mesoscale.js';
import { MolecularSurfacePresetEntry } from '../representation-presets/molecular-surface.js';
import { PolymerAndLigandPresetEntry } from '../representation-presets/polymer-and-ligand.js';
import { PolymerCartoonPresetEntry } from '../representation-presets/polymer-cartoon.js';
import { ProteinAndNucleicPresetEntry } from '../representation-presets/protein-and-nucleic.js';
import { AllModelsHierarchyPresetEntry } from '../hierarchy-presets/all-models.js';
import { CrystalContactsHierarchyPresetEntry } from '../hierarchy-presets/crystal-contacts.js';
import { DefaultHierarchyPresetEntry } from '../hierarchy-presets/default.js';
import { SupercellHierarchyPresetEntry } from '../hierarchy-presets/supercell.js';
import { UnitcellHierarchyPresetEntry } from '../hierarchy-presets/unitcell.js';

const Fixture = [
  'ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N',
  'ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 10.00           C',
  'ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 10.00           C',
  'ATOM      4  O   ALA A   1       1.300   2.400   0.000  1.00 10.00           O',
  'ATOM      5  N   GLY A   2       3.300   1.400   0.000  1.00 10.00           N',
  'ATOM      6  CA  GLY A   2       4.000   2.700   0.000  1.00 10.00           C',
  'ATOM      7  C   GLY A   2       5.500   2.500   0.000  1.00 10.00           C',
  'ATOM      8  O   GLY A   2       6.000   1.400   0.000  1.00 10.00           O',
  'TER       9      GLY A   2',
  'HETATM   10  C1  LIG B   1      10.000  10.000  10.000  1.00 10.00           C',
  'HETATM   11  O1  LIG B   1      11.200  10.000  10.000  1.00 10.00           O',
  'HETATM   12  O   HOH C   1      20.000  20.000  20.000  1.00 10.00           O',
  'HETATM   13 NA    NA D   1      30.000  30.000  30.000  1.00 10.00          NA',
  'END',
].join('\n');

async function createPlugin(entry: PluginRegistryEntry) {
  const plugin = new PluginContext({ behaviors: [], registry: [] });
  const warnings: string[] = [];
  plugin.events.log.subscribe((e) => {
    if (e.type === 'warning') warnings.push(e.message);
  });
  await plugin.init();

  // The registries preload the built-in providers; only what the entry brings may be there.
  for (const scope of [
    plugin.representation.structure,
    plugin.representation.volume,
    plugin.representation.particles,
  ]) {
    scope.registry.clear();
    scope.themes.colorThemeRegistry.clear();
    scope.themes.sizeThemeRegistry.clear();
  }
  plugin.register(entry);

  // Presets build and commit a representation tree. Nothing here renders, so the commits are skipped.
  const build = plugin.state.data.build.bind(plugin.state.data);
  plugin.state.data.build = () => {
    const builder = build();
    builder.commit = async () => {};
    return builder;
  };
  return { plugin, warnings };
}

async function loadStructure(plugin: PluginContext) {
  const data = await plugin.builders.data.rawData({ data: Fixture, label: 'fixture' });
  const trajectory = await plugin.builders.structure.parseTrajectory(data, 'pdb');
  const model = await plugin.builders.structure.createModel(trajectory);
  const structure = await plugin.builders.structure.createStructure(model);
  return { trajectory, structure };
}

/** Records the props every preset passes to `buildRepresentation`. */
function recordRepresentations(plugin: PluginContext) {
  const builder = plugin.builders.structure.representation;
  const calls: any[] = [];
  const buildRepresentation = builder.buildRepresentation.bind(builder);
  builder.buildRepresentation = ((...args: any[]) => {
    calls.push(args[2]);
    return (buildRepresentation as any)(...args);
  }) as any;
  return calls;
}

const RepresentationEntries: [string, PluginRegistryEntry, keyof typeof PresetStructureRepresentations][] = [
  ['empty', EmptyPresetEntry, 'empty'],
  ['auto', AutoPresetEntry, 'auto'],
  ['atomic-detail', AtomicDetailPresetEntry, 'atomic-detail'],
  ['polymer-cartoon', PolymerCartoonPresetEntry, 'polymer-cartoon'],
  ['polymer-and-ligand', PolymerAndLigandPresetEntry, 'polymer-and-ligand'],
  ['protein-and-nucleic', ProteinAndNucleicPresetEntry, 'protein-and-nucleic'],
  ['coarse-surface', CoarseSurfacePresetEntry, 'coarse-surface'],
  ['illustrative', IllustrativePresetEntry, 'illustrative'],
  ['molecular-surface', MolecularSurfacePresetEntry, 'molecular-surface'],
  ['auto-lod', AutoLodPresetEntry, 'auto-lod'],
  ['mesoscale', MesoscalePresetEntry, 'mesoscale'],
];

const HierarchyEntries: [PluginRegistryEntry, keyof typeof PresetTrajectoryHierarchy][] = [
  [DefaultHierarchyPresetEntry, 'default'],
  [AllModelsHierarchyPresetEntry, 'all-models'],
  [UnitcellHierarchyPresetEntry, 'unitcell'],
  [SupercellHierarchyPresetEntry, 'supercell'],
  [CrystalContactsHierarchyPresetEntry, 'crystalContacts'],
];

/** Thresholds that make the fixture small, medium, large, huge, and gigantic. */
const SizeThresholds: Partial<Structure.SizeThresholds>[] = [
  {},
  { smallResidueCount: 1 },
  { smallResidueCount: 1, mediumResidueCount: 2 },
  { smallResidueCount: 1, mediumResidueCount: 2, largeResidueCount: 2, highSymmetryUnitCount: 0 },
  { smallResidueCount: 1, mediumResidueCount: 2, largeResidueCount: 2, highSymmetryUnitCount: 1000 },
];

describe('representation preset entries', () => {
  it('list their own preset first and the built-in preset under its catalog key', () => {
    for (const [, entry, key] of RepresentationEntries) {
      expect(Object.keys(entry)).toEqual(['structure']);
      expect(Object.keys(entry.structure!).sort()).toEqual(
        // the empty preset builds nothing
        key === 'empty' ? ['presets'] : ['presets', 'representations', 'themes'],
      );
      expect(entry.structure!.presets!.representation![0]).toBe(PresetStructureRepresentations[key]);
    }
  });

  it('the auto entry includes the presets it chooses between', () => {
    const ids = AutoPresetEntry.structure!.presets!.representation!.map((p) => p.id);
    expect(ids).toEqual([
      'preset-structure-representation-auto',
      'preset-structure-representation-coarse-surface',
      'preset-structure-representation-polymer-cartoon',
      'preset-structure-representation-polymer-and-ligand',
      'preset-structure-representation-atomic-detail',
    ]);
  });

  for (const [name, entry] of RepresentationEntries) {
    it(`'${name}' registered alone has everything its preset builds`, async () => {
      const { plugin, warnings } = await createPlugin(entry);
      const calls = recordRepresentations(plugin);
      const { structure } = await loadStructure(plugin);
      const preset = plugin.builders.structure.representation.resolveProvider(name)!;
      expect(preset).toBeDefined();

      const { registry, themes } = plugin.representation.structure;
      for (const thresholds of SizeThresholds) {
        calls.length = 0;
        plugin.config.set(PluginConfig.Structure.SizeThresholds, { ...Structure.DefaultSizeThresholds, ...thresholds });
        await plugin.builders.structure.representation.applyPreset(structure, preset);

        if (name !== 'empty') expect(calls.length).toBeGreaterThan(0);
        for (const props of calls) {
          // Providers are used as they are, so each must be the one that is registered.
          expect(registry.has(props.type)).toBe(true);
          const type = typeof props.type === 'string' ? registry.get(props.type) : props.type;
          expect(themes.colorThemeRegistry.has(props.color ?? type.defaultColorTheme.name)).toBe(true);
          expect(themes.sizeThemeRegistry.has(props.size ?? type.defaultSizeTheme.name)).toBe(true);
        }
      }
      expect(warnings.filter((w) => /not registered/.test(w))).toEqual([]);
      plugin.dispose();
    });
  }

  it('the preset modules do not import catalogs', () => {
    const dir = path.resolve(__dirname, '../representation-presets');
    let count = 0;
    for (const file of fs.readdirSync(dir)) {
      if (file === 'catalog.ts') continue;
      const imports = fs
        .readFileSync(path.join(dir, file), 'utf8')
        .split('\n')
        .filter((l) => l.startsWith('import '));
      expect(imports.filter((l) => /catalog/.test(l))).toEqual([]);
      count++;
    }
    expect(count).toBe(RepresentationEntries.length + 1); // + types
  });

  it('DefaultPresets lists the preset of each entry in the 5.x order', () => {
    expect(DefaultPresets.structure!.presets!.representation).toEqual(Object.values(PresetStructureRepresentations));
    expect(DefaultPresets.structure!.presets!.hierarchy).toEqual(Object.values(PresetTrajectoryHierarchy));
    expect(DefaultPresets.structure!.presets!.representation!.length).toBe(11);
    expect(DefaultPresets.structure!.presets!.hierarchy!.length).toBe(5);
  });
});

describe('hierarchy preset entries', () => {
  it('list their own preset first', () => {
    for (const [entry, key] of HierarchyEntries) {
      expect(entry.structure!.presets!.hierarchy![0]).toBe(PresetTrajectoryHierarchy[key]);
    }
  });

  it('do not include representation presets', () => {
    for (const [entry] of HierarchyEntries) {
      expect(entry.structure!.presets!.representation).toBeUndefined();
      expect(entry.structure!.representations).toBeUndefined();
    }
  });

  it('the hierarchy preset modules do not import the representation preset catalog or the auto preset', () => {
    const dir = path.resolve(__dirname, '../hierarchy-presets');
    for (const file of fs.readdirSync(dir)) {
      if (file === 'catalog.ts') continue;
      const imports = fs
        .readFileSync(path.join(dir, file), 'utf8')
        .split('\n')
        .filter((l) => l.startsWith('import ') && !l.startsWith('import type'));
      expect(imports.filter((l) => /representation-presets\/(auto|catalog)/.test(l))).toEqual([]);
    }
  });
});

describe('hierarchy presets', () => {
  const ConfiguredId = 'test-configured-representation-preset';
  const ParamId = 'test-param-representation-preset';

  async function setup(entry: PluginRegistryEntry) {
    const { plugin, warnings } = await createPlugin(entry);
    const applied: { id: string; params: any }[] = [];
    const preset = (id: string) => ({
      id,
      display: { name: id },
      apply: async (_: any, params: any) => {
        applied.push({ id, params });
        return {};
      },
    });
    plugin.register({
      structure: { presets: { representation: [preset(ConfiguredId), preset(ParamId)] } },
    });
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, ConfiguredId);
    const { trajectory } = await loadStructure(plugin);
    return { plugin, trajectory, applied, warnings };
  }

  for (const [entry, key] of HierarchyEntries) {
    it(`'${key}' applies the configured representation preset`, async () => {
      const { plugin, trajectory, applied } = await setup(entry);
      await plugin.builders.structure.hierarchy.applyPreset(trajectory, key as any);
      expect(applied.map((a) => a.id)).toEqual([ConfiguredId]);
      plugin.dispose();
    });

    it(`'${key}' prefers the representation preset it is given`, async () => {
      const { plugin, trajectory, applied } = await setup(entry);
      await plugin.builders.structure.hierarchy.applyPreset(trajectory, key as any, {
        representationPreset: ParamId as any,
      });
      expect(applied.map((a) => a.id)).toEqual([ParamId]);
      plugin.dispose();
    });
  }

  it('fails with the preset error when the configured preset is not registered', async () => {
    const { plugin, trajectory } = await setup(DefaultHierarchyPresetEntry);
    plugin.config.set(PluginConfig.Structure.DefaultRepresentationPreset, 'not-registered');
    await expect(plugin.builders.structure.hierarchy.applyPreset(trajectory, 'default')).rejects.toThrow(
      "Preset 'not-registered' is not registered in this plugin",
    );
    plugin.dispose();
  });

  it('has no built-in default for the representation preset param', async () => {
    const { plugin, trajectory } = await setup(DefaultHierarchyPresetEntry);
    const params = PresetTrajectoryHierarchy.default.params!(trajectory.cell?.obj, plugin) as any;
    expect(params.representationPreset.defaultValue).toBe('');
    plugin.dispose();
  });
});
