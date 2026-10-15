/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { PluginBehavior } from '@molstar/plugin/behavior/behavior';
import { PluginSpec } from '@molstar/plugin/spec';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';
import { AssemblySymmetryClusterColorThemeProvider } from './color.js';
import { AssemblySymmetryApi } from './operations.js';
import { AssemblySymmetryProvider } from './prop.js';

declare module '@molstar/prism-viewer/api' {
  // The `assembly-symmetry` area (parity.md §6): RCSB PDB's assembly symmetry as operations.
  interface ExtensionApis {
    'assembly-symmetry': AssemblySymmetryApi;
  }
}

// Assembly symmetry in the Prism viewer (plan 5.2, README decision 27): commands first. The symmetry RCSB PDB annotates
// for an assembly (`get`), its axes (`show`, `hide`) and the coloring by its clusters (`color`); the plugin gets the
// property and the cluster color theme. v5's preset and state action are left out: they fetch from RCSB PDB outside
// the operations, past the agent's approval (agent.md §4). Its view is a row under each structure that is a PDB
// entry's assembly (`assembly-symmetry.row.symmetry`, parity.md §6; viewer-spec.md §3, `rows`), in a chunk of its own.

/**
 * The cluster theme, which never fetches by itself: it colors by what `show` or `color` (or the structure's
 * custom-properties node, on a session) fetched, so `representation.color` with it reaches no network.
 */
export const ClusterColorTheme: typeof AssemblySymmetryClusterColorThemeProvider = {
  ...AssemblySymmetryClusterColorThemeProvider,
  ensureCustomProperties: {
    attach: async (ctx, data) => {
      const s = data.structure?.root;
      if (s && AssemblySymmetryProvider.get(s).value) {
        await AssemblySymmetryClusterColorThemeProvider.ensureCustomProperties!.attach(ctx, data);
      }
    },
    detach: (data) => AssemblySymmetryClusterColorThemeProvider.ensureCustomProperties!.detach(data),
  },
};

/**
 * The property registered for the structures' custom-properties node, without v5's React controls (`behavior.ts`
 * sets those): the node then holds the chosen symmetry.
 */
export const AssemblySymmetryProperty = PluginBehavior.create<{ autoAttach: boolean }>({
  name: 'prism-assembly-symmetry-prop',
  category: 'custom-props',
  display: {
    name: 'Assembly Symmetry',
    description: 'Assembly Symmetry data provided by RCSB PDB (calculated with BioJava) or by PDBe.',
  },
  ctor: class extends PluginBehavior.Handler<{ autoAttach: boolean }> {
    register(): void {
      this.ctx.customStructureProperties.register(AssemblySymmetryProvider, this.params.autoAttach);
    }
    unregister() {
      this.ctx.customStructureProperties.unregister(AssemblySymmetryProvider.descriptor.name);
    }
  },
  params: () => ({ autoAttach: PD.Boolean(false) }),
});

export const AssemblySymmetrySkill: ViewerSkill = {
  id: 'assembly-symmetry',
  title: 'Assembly symmetry',
  summary: "a PDB assembly's symmetry from RCSB PDB: its point group, axes, coloring by cluster",
  feature: 'assembly-symmetry',
  guidance: [
    "Skill assembly-symmetry: RCSB PDB annotates the symmetry of a PDB entry's assemblies (global, pseudo, local).",
    'get_assembly_symmetry lists them with their number, point group (C2, D3, …), stoichiometry and clusters;',
    'show_symmetry_axes draws the axes of one (0: the first that is not C1), hide_symmetry_axes removes them;',
    'color_by_symmetry colors the chains by the symmetry cluster they belong to.',
    'Each asks RCSB PDB the first time; a structure not loaded from the PDB has none.',
  ].join('\n'),
  tools: [
    {
      id: 'get_assembly_symmetry',
      operation: 'assembly-symmetry.get',
      description: "Lists the symmetries RCSB PDB has for a structure's assembly.",
    },
    {
      id: 'show_symmetry_axes',
      operation: 'assembly-symmetry.show',
      description: "Shows the symmetry axes of a structure's assembly.",
    },
    { id: 'hide_symmetry_axes', operation: 'assembly-symmetry.hide', description: 'Removes the symmetry axes.' },
    {
      id: 'color_by_symmetry',
      operation: 'assembly-symmetry.color',
      description: 'Colors a structure by symmetry cluster (chains related by the symmetry share a color).',
    },
  ],
};

export const AssemblySymmetryExtension: ViewerExtension = {
  id: 'assembly-symmetry',
  title: 'Assembly symmetry',
  plugin: {
    registry: [{ structure: { themes: { color: [ClusterColorTheme] } } }],
    behaviors: [PluginSpec.Behavior(AssemblySymmetryProperty)],
  },
  api: { 'assembly-symmetry': (ctx) => new AssemblySymmetryApi(ctx) },
  commands: (api) => {
    const a = api['assembly-symmetry'];
    if (!a) return [];
    const category = 'Structure';
    return [
      operationCommand(api, a.show, {
        category,
        icon: 'orbit',
        alias: ['symmetry', 'symmetry axes', 'point group', 'rcsb symmetry'],
      }),
      operationCommand(api, a.hide, { category, icon: 'eye-off', alias: ['hide symmetry'] }),
      operationCommand(api, a.color, {
        category,
        icon: 'paintbrush',
        alias: ['symmetry cluster', 'color by symmetry'],
      }),
      operationCommand(api, a.get, { category, icon: 'info', alias: ['stoichiometry', 'what symmetry'] }),
    ];
  },
  skills: [AssemblySymmetrySkill],
  ui: async (ctx) => (await import('./prism-viewer-ui.js')).assemblySymmetryUI(ctx),
};
