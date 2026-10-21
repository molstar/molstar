/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { CommandDefinition } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';
import { operationCommand } from './operations.js';

/** A param kept out of the sentence: set by the command's own form in the command menu (`ui/command-forms.ts`), if any. */
const hidden = (params: PD.Params, ...names: string[]): PD.Params => {
  const out = { ...params };
  for (const k of names) if (out[k]) out[k] = { ...out[k], isHidden: true };
  return out;
};

/** The representation commands (api.md §9, §2 `representation`): the look of a target, and what a component shows. */
export function representationCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  const category = 'Style';
  return [
    operationCommand(api, api.representation.color, {
      id: 'representation.color',
      title: 'Color…',
      category,
      icon: 'paintbrush',
      alias: ['color', 'paint', 'colour', 'recolor'],
      action: 'Color',
      values: { theme: '' },
    }),
    operationCommand(api, api.representation.color, {
      id: 'representation.color-theme',
      title: 'Color by theme…',
      category,
      icon: 'palette',
      alias: ['color theme', 'color by', 'theme', 'color scheme'],
      action: 'Apply',
      // A sentence: "ligand in 1TQN by Chain Id".
      order: ['target', 'structure', 'theme'],
      fields: { structure: { before: 'in' }, theme: { before: 'by' } },
      // A theme is asked for (no "None"); the color is the uniform theme's, set by the form's color picker.
      form: (params) => {
        const out = hidden(params, 'color');
        const theme = out.theme as PD.Select<string>;
        if (theme) out.theme = { ...theme, options: theme.options.filter(([v]) => v !== '') };
        return out;
      },
    }),
    operationCommand(api, api.representation.transparency, {
      category,
      icon: 'eye-off',
      alias: ['transparent', 'opacity', 'alpha', 'ghost', 'make transparent'],
      action: 'Set',
    }),
    operationCommand(api, api.representation.set, {
      category,
      icon: 'layers',
      title: 'Show as…',
      alias: [
        'show as',
        'representation',
        'cartoon',
        'ball and stick',
        'spacefill',
        'surface',
        'change representation',
      ],
      action: 'Show',
      // A sentence: "Polymer in 1TQN as Cartoon".
      order: ['component', 'structure', 'type'],
      fields: { structure: { before: 'in' }, type: { before: 'as' } },
      // The type's options are the form's (size, quality, …); typed, the type's defaults.
      form: (params) => hidden(params, 'params'),
    }),
    operationCommand(api, api.representation.remove, {
      category,
      icon: 'trash-2',
      title: 'Remove a representation…',
      alias: ['remove representation', 'hide representation', 'delete representation'],
      action: 'Remove',
      // By component and type here; one by state ref is the tree row's.
      values: { node: '' },
    }),
    operationCommand(api, api.representation.visibility, {
      category,
      icon: 'eye',
      alias: ['show representation', 'toggle representation'],
      action: 'Set',
    }),
    operationCommand(api, api.representation.size, {
      category,
      icon: 'maximize-2',
      alias: ['size', 'thickness', 'radius', 'scale', 'size theme'],
      action: 'Set',
    }),
    operationCommand(api, api.representation.emissive, {
      category,
      icon: 'sun',
      alias: ['emissive', 'glow', 'emission'],
      action: 'Set',
    }),
    operationCommand(api, api.representation.material, {
      category,
      icon: 'sparkles',
      alias: ['material', 'metalness', 'roughness', 'shiny', 'matte'],
      action: 'Set',
    }),
    operationCommand(api, api.representation.clipping, {
      category,
      icon: 'crop',
      alias: ['clipping', 'exclude from clipping', 'clip', 'clip objects'],
      action: 'Set',
    }),
    operationCommand(api, api.representation.reset, {
      category,
      icon: 'rotate-ccw',
      alias: ['reset color', 'reset transparency', 'clear color', 'clear style', 'reset style', 'reset emissive'],
      action: 'Reset',
    }),
    operationCommand(api, api.representation.layerUpdate, {
      category,
      icon: 'pencil',
      title: 'Change a style layer…',
      alias: ['edit color', 'change color', 'edit transparency', 'edit glow', 'edit layer', 'edit overpaint'],
      action: 'Change',
    }),
    operationCommand(api, api.representation.layerRemove, {
      category,
      icon: 'trash-2',
      title: 'Remove a style layer…',
      alias: ['remove color', 'delete color', 'remove overpaint', 'remove transparency', 'remove glow', 'remove layer'],
      action: 'Remove',
    }),
  ];
}
