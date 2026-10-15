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
    operationCommand(api, api.representation.reset, {
      category,
      icon: 'rotate-ccw',
      alias: ['reset color', 'reset transparency', 'clear color', 'clear style', 'reset style', 'reset emissive'],
      action: 'Reset',
    }),
  ];
}
