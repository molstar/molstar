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

/** The component commands (api.md §9, §2 `component`): one per operation, the reads and `modify` aside. */
export function componentCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  const category = 'Components';
  return [
    operationCommand(api, api.component.add, {
      category,
      icon: 'plus',
      alias: ['add component', 'new component', 'make component', 'show ligand', 'add representation'],
      action: 'Add',
    }),
    operationCommand(api, api.component.modify, {
      category,
      icon: 'layers',
      alias: ['modify component', 'edit component', 'add to component', 'remove from component'],
      action: 'Change',
    }),
    operationCommand(api, api.component.rename, {
      category,
      icon: 'pencil',
      alias: ['rename component'],
      action: 'Rename',
    }),
    operationCommand(api, api.component.visibility, {
      category,
      icon: 'eye',
      alias: ['show component', 'hide component', 'toggle component'],
      action: 'Set',
    }),
    operationCommand(api, api.component.remove, {
      category,
      icon: 'trash-2',
      alias: ['delete component', 'remove component'],
      action: 'Remove',
    }),
    operationCommand(api, api.look.set, {
      category: 'Style',
      icon: 'sliders-horizontal',
      title: 'Set how representations are drawn…',
      alias: ['scene defaults', 'component options', 'hydrogens', 'quality', 'material', 'ignore light', 'look'],
      action: 'Set',
    }),
    operationCommand(api, api.look.reset, {
      category: 'Style',
      icon: 'rotate-ccw',
      title: 'Use the scene defaults again…',
      alias: ['reset look', 'clear overrides', 'scene defaults'],
      action: 'Reset',
    }),
    operationCommand(api, api.look.get, {
      category: 'Style',
      icon: 'list',
      alias: ['overrides', 'look', 'scene defaults'],
    }),
    operationCommand(api, api.contents.visibility, {
      category,
      icon: 'eye',
      title: 'Show or hide an entity or chain…',
      alias: ['show entity', 'hide entity', 'show chain', 'hide chain', 'toggle chain'],
      action: 'Set',
    }),
    operationCommand(api, api.component.list, {
      category,
      icon: 'list',
      alias: ['components', 'list components'],
    }),
  ];
}
