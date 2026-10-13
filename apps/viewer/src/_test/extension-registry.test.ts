/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ANVILMembraneOrientation } from '@molstar/anvil-extension/behavior';
import { AssemblySymmetry } from '@molstar/assembly-symmetry-extension';
import { InitAssemblySymmetry3D } from '@molstar/assembly-symmetry-extension/behavior';
import { DnatcoNtCs } from '@molstar/dnatco-extension';
import { G3DFormat, LoadG3D } from '@molstar/g3d-extension/format';
import { KinemageExtension } from '@molstar/kinemage-extension/behavior';
import { MAQualityAssessment } from '@molstar/model-archive-extension/quality-assessment/behavior';
import { PDBeStructureQualityReport } from '@molstar/pdbe-extension';
import { RCSBValidationReport } from '@molstar/rcsb-extension';
import { SbNcbrPartialCharges, SbNcbrTunnels } from '@molstar/sb-ncbr-extension';
import { DownloadTunnels } from '@molstar/sb-ncbr-extension/tunnels/actions';
import { wwPDBChemicalComponentDictionary } from '@molstar/wwpdb-extension/ccd/behavior';
import type { StateAction, StateTransformer } from '@molstar/core/state';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { PluginContext } from '@molstar/plugin/context';
import { DefaultPluginSpec } from '@molstar/plugin/default-spec';
import { createViewerRegistry } from '@molstar/viewer/registry';

/**
 * The Viewer extensions whose behaviors register their providers (representations, themes, presets, selection
 * queries, loci labels, formats, drag-and-drop handlers and actions) through one `plugin.register(entry)` call.
 */
const ConvertedBehaviors = [
  DnatcoNtCs,
  ANVILMembraneOrientation,
  AssemblySymmetry,
  RCSBValidationReport,
  PDBeStructureQualityReport,
  MAQualityAssessment,
  SbNcbrPartialCharges,
  SbNcbrTunnels,
  KinemageExtension,
  G3DFormat,
  wwPDBChemicalComponentDictionary,
];

async function createPlugin(behaviors: StateTransformer[]) {
  const defaults = DefaultPluginSpec();
  const spec = {
    ...defaults,
    registry: createViewerRegistry(),
    behaviors: [...defaults.behaviors, ...behaviors.map((b) => ({ transformer: b, defaultParams: {} }))],
  };
  const plugin = new PluginContext(spec as any);
  await plugin.init();
  return plugin;
}

/** What the registries, managers, and custom property registries hold, in registration order. */
function snapshot(plugin: PluginContext) {
  const names = (list: readonly { name: string }[]) => list.map((e) => e.name);
  const { structure, volume, particles } = plugin.representation;
  const scope = (s: typeof structure) => ({
    representations: names(s.registry.list),
    color: names(s.themes.colorThemeRegistry.list),
    size: names(s.themes.sizeThemeRegistry.list),
  });
  const actions = (plugin.state.data.actions as any).actions as Map<string, { action: StateAction; count: number }>;
  return {
    structure: scope(structure),
    volume: scope(volume as typeof structure),
    particles: scope(particles as typeof structure),
    formats: names(plugin.dataFormats.list),
    hierarchyPresets: plugin.builders.structure.hierarchy.providers.map((p) => p.id),
    representationPresets: plugin.builders.structure.representation.providers.map((p) => p.id),
    queries: plugin.query.structure.registry.list.map((q) => `${q.category}: ${q.label}`),
    lociLabels: plugin.managers.lociLabels.providers.length,
    markdown: names(plugin.managers.markdownExtensions.extension),
    dragAndDrop: plugin.managers.dragAndDrop.list().map((e) => e.name),
    animations: names(plugin.managers.animation.animations),
    actions: [...actions.values()].map((e) => `${e.action.definition.display.name} x${e.count}`),
    customProperties: {
      model: [...plugin.customModelProperties.providers.keys()],
      structure: [...plugin.customStructureProperties.providers.keys()],
      volume: [...plugin.customVolumeProperties.providers.keys()],
    },
  };
}

/** The snapshot with every list sorted, to compare registries whose provider order changes when re-registered. */
function sorted<T>(value: T): T {
  if (Array.isArray(value)) return [...value].sort() as T;
  if (value && typeof value === 'object') {
    return Object.fromEntries(Object.entries(value).map(([k, v]) => [k, sorted(v)])) as T;
  }
  return value;
}

/** The order in which `PluginContext.init` creates behaviors: custom properties first, then by category. */
function initializationOrder(behaviors: StateTransformer[]) {
  const categories = Object.keys(PluginBehavior.Categories);
  const rank = (b: StateTransformer) => {
    const category = PluginBehavior.getCategoryId(b);
    return category === 'custom-props' ? -1 : categories.indexOf(category);
  };
  return behaviors
    .map((b, i) => ({ b, i }))
    .sort((x, y) => rank(x.b) - rank(y.b) || x.i - y.i)
    .map(({ b }) => b);
}

function disable(plugin: PluginContext, behavior: StateTransformer) {
  return plugin.runTask(plugin.state.behaviors.updateTree(plugin.state.behaviors.build().delete(behavior.id)));
}

function enable(plugin: PluginContext, behavior: StateTransformer) {
  return plugin.state.updateBehavior(behavior, () => {});
}

function isEnabled(plugin: PluginContext, behavior: StateTransformer) {
  return plugin.state.behaviors.tree.transforms.has(behavior.id);
}

function actionCount(plugin: PluginContext, action: StateAction) {
  const actions = (plugin.state.data.actions as any).actions as Map<string, { count: number }>;
  return actions.get(action.id)?.count ?? 0;
}

describe('extension behaviors registering through plugin.register', () => {
  let withExtensions: PluginContext;
  let initial: ReturnType<typeof snapshot>;

  beforeAll(async () => {
    withExtensions = await createPlugin(ConvertedBehaviors);
    initial = snapshot(withExtensions);
  });

  afterAll(() => withExtensions.dispose());

  it('registers the extension providers on top of the Viewer registry', async () => {
    const base = snapshot(await createPlugin([]));
    expect(initial.structure.representations.length).toBeGreaterThan(base.structure.representations.length);
    expect(initial.structure.color.length).toBeGreaterThan(base.structure.color.length);
    expect(initial.representationPresets.length).toBeGreaterThan(base.representationPresets.length);
    expect(initial.hierarchyPresets).toContain('preset-trajectory-ccd');
    expect(initial.formats).toContain('KIN');
    expect(initial.dragAndDrop).toContain('kin');
    expect(initial.lociLabels).toBeGreaterThan(base.lociLabels);
  });

  for (const behavior of ConvertedBehaviors) {
    it(`leaves the registries as before when ${behavior.id} is toggled off and on`, async () => {
      expect(sorted(snapshot(withExtensions))).toEqual(sorted(initial));

      await disable(withExtensions, behavior);
      expect(isEnabled(withExtensions, behavior)).toBe(false);
      const off = snapshot(withExtensions);
      // the behavior owned something: at least one registry shrank
      expect(sorted(off)).not.toEqual(sorted(initial));
      const lost = JSON.stringify(sorted(initial)).length - JSON.stringify(sorted(off)).length;
      expect(lost).toBeGreaterThan(0);

      await enable(withExtensions, behavior);
      expect(isEnabled(withExtensions, behavior)).toBe(true);
      expect(sorted(snapshot(withExtensions))).toEqual(sorted(initial));

      // toggling again does not accumulate registrations
      await disable(withExtensions, behavior);
      await enable(withExtensions, behavior);
      expect(sorted(snapshot(withExtensions))).toEqual(sorted(initial));
    });
  }

  it('removes everything the behaviors registered, and restores it in the original order', async () => {
    const without = await createPlugin([]);
    const baseline = snapshot(without);
    without.dispose();

    const plugin = await createPlugin(ConvertedBehaviors);
    expect(snapshot(plugin)).toEqual(initial);

    for (const b of ConvertedBehaviors) await disable(plugin, b);
    const off = snapshot(plugin);
    // behaviors that stay (and the global custom property registries of the extensions) aside, the registries match
    expect({ ...off, customProperties: undefined }).toEqual({ ...baseline, customProperties: undefined });

    for (const b of initializationOrder(ConvertedBehaviors)) await enable(plugin, b);
    expect(snapshot(plugin)).toEqual(initial);
    plugin.dispose();
  });

  describe('actions are counted', () => {
    const cases = [
      ['assembly symmetry', AssemblySymmetry, InitAssemblySymmetry3D],
      ['SB-NCBR tunnels', SbNcbrTunnels, DownloadTunnels],
      ['g3d', G3DFormat, LoadG3D],
    ] as const;

    for (const [name, behavior, action] of cases) {
      it(`${name}: removing the behavior does not remove an action the spec also lists`, async () => {
        const plugin = await createPlugin([behavior]);
        expect(actionCount(plugin, action)).toBe(1);

        // the app lists the same action as well, as a registry entry
        const undoSpecEntry = plugin.register({ actions: [action] });
        expect(actionCount(plugin, action)).toBe(2);

        await disable(plugin, behavior);
        expect(actionCount(plugin, action)).toBe(1);
        const byType = (plugin.state.data.actions as any).fromTypeIndex as Map<unknown, StateAction[]>;
        expect(byType.get(action.definition.from[0].type)).toContain(action);

        await enable(plugin, behavior);
        expect(actionCount(plugin, action)).toBe(2);

        undoSpecEntry();
        expect(actionCount(plugin, action)).toBe(1);
        await disable(plugin, behavior);
        expect(actionCount(plugin, action)).toBe(0);
        plugin.dispose();
      });
    }
  });
});
