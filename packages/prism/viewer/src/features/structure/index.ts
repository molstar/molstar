/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '@molstar/prism-ui/core';
import { ComponentApi } from '../../api/component.js';
import { ContentsApi } from '../../api/contents.js';
import { LookApi } from '../../api/look.js';
import type { ViewerCommandContext } from '../../api/managers.js';
import { RepresentationApi } from '../../api/representation.js';
import { StructureApi } from '../../api/structure.js';
import { StylesApi } from '../../api/styles.js';
import { componentCommands } from '../../commands/component.js';
import { operationCommand } from '../../commands/operations.js';
import { representationCommands } from '../../commands/representation.js';
import type { ViewerApi } from '../../api/area.js';
import type { ViewerExtension } from '../../spec.js';

// The `structure` feature (viewer-spec.md §3.1): the scene's API, without its panel: `structure` (presets, models,
// assemblies, symmetry, transforms), `contents`, `component`, `look`, `representation`, `styles`, and their commands. The
// Structure panel is the `structure-panel` feature, on top of this one.

export const StructureFeature: ViewerExtension = {
  id: 'structure',
  title: 'Structure',
  builtIn: true,
  api: {
    structure: (c) => new StructureApi(c),
    contents: (c) => new ContentsApi(c),
    component: (c) => new ComponentApi(c),
    look: (c) => new LookApi(c),
    representation: (c) => new RepresentationApi(c),
    styles: (c) => new StylesApi(c),
  },
  commands: (api) => structureFeatureCommands(api),
};

function structureFeatureCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  return [
    // The default palette (the representation area's, this feature's): with the representations it colors.
    operationCommand(api, api.representation.defaultPalette, {
      title: 'Use colorblind-safe colors',
      category: 'Canvas',
      icon: 'eye',
      alias: ['colorblind', 'color blind', 'okabe ito', 'palette', 'accessible colors'],
    }),
    operationCommand(api, api.structure.get, {
      category: 'Data',
      icon: 'list',
      alias: ['what is loaded', 'scene', 'outline'],
    }),
    operationCommand(api, api.structure.preset, {
      category: 'Data',
      icon: 'sparkles',
      title: 'Apply a representation preset…',
      alias: ['preset', 'representation preset', 'polymer cartoon', 'illustrative', 'auto'],
      action: 'Apply',
    }),
    operationCommand(api, api.structure.assembly, {
      category: 'Data',
      icon: 'box',
      alias: ['assembly', 'biological assembly', 'unit cell', 'deposited coordinates'],
      action: 'Show',
    }),
    operationCommand(api, api.structure.visibility, {
      category: 'Data',
      icon: 'eye',
      alias: ['show structure', 'hide structure', 'toggle structure'],
      action: 'Set',
    }),
    operationCommand(api, api.structure.remove, {
      category: 'Data',
      icon: 'trash-2',
      alias: ['delete structure', 'remove structure', 'close structure'],
      action: 'Remove',
    }),
    operationCommand(api, api.styles.apply, {
      category: 'Style',
      icon: 'swatch-book',
      alias: ['quick style', 'look', 'style', 'default style', 'illustrative style'],
      action: 'Apply',
    }),
    operationCommand(api, api.structure.bonds, {
      category: 'Data',
      icon: 'ball-and-stick',
      alias: ['bonds', 'bond provider', 'bond orders', 'no computed bonds', 'connectivity'],
      action: 'Choose',
    }),
    operationCommand(api, api.structure.symmetry, {
      category: 'Data',
      icon: 'hexagon',
      alias: [
        'symmetry',
        'symmetry mates',
        'crystal',
        'unit cell',
        'lattice',
        'periodic images',
        'neighbouring boxes',
        'tile box',
        'supercell',
      ],
    }),
    operationCommand(api, api.structure.unitCell, {
      category: 'Data',
      icon: 'box',
      alias: ['unit cell', 'crystal cell', 'cell box', 'lattice box'],
      action: 'Show',
    }),
    operationCommand(api, api.structure.transform, {
      category: 'Data',
      icon: 'move-3d',
      alias: ['transform', 'move', 'rotate', 'translate', 'conformation transform', 'matrix'],
      action: 'Move',
    }),
    operationCommand(api, api.structure.transformReset, {
      category: 'Data',
      icon: 'rotate-ccw',
      alias: ['reset transform', 'reset position', 'untransform'],
      action: 'Reset',
    }),
    ...componentCommands(api),
    ...representationCommands(api),
  ];
}
