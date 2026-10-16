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

/**
 * The selection commands (api.md §9): select (the form has the mode), the selection's edits, the focus and the
 * camera on a target. Each runs its operation as the person; none is a history step (sequence.md §7.3).
 */
export function selectionCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  return [
    operationCommand(api, api.select.set, {
      category: 'Selection',
      icon: 'mouse-pointer-2',
      title: 'Select…',
      alias: ['select', 'pick', 'add to selection', 'remove from selection'],
      action: 'Select',
      // In the Select group: its mode starts from the group's (select.ts).
      group: 'select',
    }),
    operationCommand(api, api.select.all, {
      category: 'Selection',
      icon: 'scan',
      alias: ['select everything'],
      values: { structure: '' },
    }),
    operationCommand(api, api.select.clear, {
      category: 'Selection',
      icon: 'x',
      alias: ['deselect', 'unselect', 'clear selection'],
    }),
    operationCommand(api, api.select.invert, {
      category: 'Selection',
      icon: 'squares-exclude',
      alias: ['inverse', 'everything else'],
      values: { structure: '' },
    }),
    operationCommand(api, api.select.around, {
      category: 'Selection',
      icon: 'squares-unite',
      title: 'Select the surroundings…',
      alias: ['around', 'within', 'neighbours', 'nearby'],
      action: 'Select',
    }),
    operationCommand(api, api.select.expand, {
      category: 'Selection',
      icon: 'maximize-2',
      title: 'Expand the selection…',
      alias: ['expand', 'whole residues', 'whole chains'],
      action: 'Expand',
    }),
    operationCommand(api, api.select.mode, {
      category: 'Selection',
      icon: 'mouse-pointer-click',
      alias: ['selection mode', 'click to select', 'pick mode'],
      action: 'Set',
    }),
    operationCommand(api, api.select.granularity, {
      category: 'Selection',
      icon: 'pick-residue',
      alias: ['picking level', 'granularity', 'pick residues', 'pick atoms', 'pick chains'],
      action: 'Set',
    }),
    operationCommand(api, api.focus.set, {
      category: 'Selection',
      icon: 'crosshair',
      title: 'Focus…',
      alias: ['focus', 'look at', 'zoom to'],
      action: 'Focus',
    }),
    operationCommand(api, api.focus.clear, {
      category: 'Selection',
      icon: 'x',
      alias: ['unfocus', 'clear focus'],
    }),
    operationCommand(api, api.camera.frame, {
      category: 'View',
      icon: 'locate',
      title: 'Frame…',
      alias: ['frame', 'bring into view', 'center on'],
      action: 'Frame',
      values: { animate: true },
    }),
  ];
}
