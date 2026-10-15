/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../api/area.js';
import type { ViewerCommandContext } from '../api/managers.js';
import { operationCommand } from './operations.js';

/** The representation commands (operations.md §3, §2 `representation`): the look of a target, and what a component shows. */
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
      values: { color: '' },
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
      values: { params: {} },
    }),
    operationCommand(api, api.representation.remove, {
      category,
      icon: 'trash-2',
      title: 'Remove a representation…',
      alias: ['remove representation', 'hide representation', 'delete representation'],
      action: 'Remove',
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
