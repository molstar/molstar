/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginSpec } from '@molstar/plugin/spec';
import { operationCommand } from '@molstar/prism-viewer/commands/operations';
import type { ViewerExtension, ViewerSkill } from '@molstar/prism-viewer/spec';
import { ANVILMembraneOrientation } from './behavior.js';
import { MembraneApi } from './operations.js';

declare module '@molstar/prism-viewer/api' {
  // The `membrane` area (parity.md §6): ANVIL's membrane orientation as operations.
  interface ExtensionApis {
    membrane: MembraneApi;
  }
}

// ANVIL's membrane orientation in the Prism viewer (parity.md §6, README decision 27): designed anew as commands. Show,
// hide and fit the membrane, set how it looks, say where it is; the plugin gets v5's behavior (the property, the
// membrane representation, its preset and the "Residues Embedded in Membrane" selection query), which brings no UI. Its
// view is a row in the Other group of each structure ANVIL applies to (`anvil.row.membrane`, parity.md §6;
// viewer-spec.md §3, `rows`), in a chunk of its own.

export const MembraneSkill: ViewerSkill = {
  id: 'membrane',
  title: 'Membrane',
  summary: 'place a membrane around a membrane protein (ANVIL), show, fit and style it',
  feature: 'anvil',
  guidance: [
    "Skill membrane: ANVIL places a lipid bilayer around a protein from its residues' hydrophobicity and exposure,",
    'computed in the viewer (it can take a moment for a large structure). show_membrane draws its two planes,',
    'hide_membrane removes them, get_membrane says its thickness, normal and how many residues are in it.',
    'fit_membrane places it again with ANVIL settings (minThickness, maxThickness, adjust: 14 for helices, 5 for',
    "beta barrels; tmdetDefinition for porins); style_membrane changes the planes' color, opacity, radius, rims.",
  ].join('\n'),
  tools: [
    { id: 'get_membrane', operation: 'membrane.get', description: 'Where the membrane of a protein is.' },
    { id: 'show_membrane', operation: 'membrane.show', description: 'Places and shows the membrane of a protein.' },
    { id: 'hide_membrane', operation: 'membrane.hide', description: 'Removes the membrane.' },
    { id: 'fit_membrane', operation: 'membrane.fit', description: 'Places the membrane again with ANVIL settings.' },
    { id: 'style_membrane', operation: 'membrane.options', description: 'Changes how the membrane looks.' },
  ],
};

export const AnvilExtension: ViewerExtension = {
  id: 'anvil',
  title: 'Membrane orientation',
  plugin: { registry: [], behaviors: [PluginSpec.Behavior(ANVILMembraneOrientation)] },
  api: { membrane: (ctx) => new MembraneApi(ctx) },
  commands: (api) => {
    const m = api.membrane;
    if (!m) return [];
    const category = 'Structure';
    return [
      operationCommand(api, m.show, {
        category,
        icon: 'layers',
        alias: ['membrane', 'anvil', 'lipid bilayer', 'membrane orientation'],
      }),
      operationCommand(api, m.hide, { category, icon: 'eye-off', alias: ['hide membrane'] }),
      operationCommand(api, m.fit, { category, icon: 'scan', alias: ['anvil', 'membrane thickness', 'refit'] }),
      operationCommand(api, m.options, {
        category,
        icon: 'sliders-horizontal',
        alias: ['membrane color', 'membrane opacity', 'membrane style'],
        action: 'Apply',
      }),
      operationCommand(api, m.get, { category, icon: 'info', alias: ['membrane thickness', 'transmembrane'] }),
    ];
  },
  skills: [MembraneSkill],
  ui: async (ctx) => (await import('./prism-viewer-ui.js')).membraneUI(ctx),
};
