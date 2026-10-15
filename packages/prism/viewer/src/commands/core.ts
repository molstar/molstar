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

// The core's commands (operations.md §3), over the API bound to the person; the features bring theirs.

export function coreCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  return [
    operationCommand(api, api.load.files, {
      category: 'Data',
      icon: 'folder-open',
      title: 'Open files…',
      alias: ['open', 'load file', 'import'],
      action: 'Open',
    }),
    operationCommand(api, api.session.open, {
      category: 'Data',
      icon: 'folder-open',
      title: 'Open a session…',
      alias: ['open session', 'molx', 'molj'],
      action: 'Open',
    }),
    operationCommand(api, api.session.clear, {
      category: 'Data',
      icon: 'trash-2',
      alias: ['clear', 'reset', 'remove everything', 'new session'],
    }),
    operationCommand(api, api.structure.get, {
      category: 'Data',
      icon: 'list',
      alias: ['what is loaded', 'scene', 'outline'],
    }),
    operationCommand(api, api.canvas.background, {
      category: 'Canvas',
      icon: 'palette',
      alias: ['background', 'bg'],
      action: 'Set',
    }),
    operationCommand(api, api.camera.reset, {
      category: 'View',
      icon: 'focus',
      alias: ['reset view', 'fit', 'recenter'],
      values: { animate: true },
    }),
    operationCommand(api, api.history.undo, {
      category: 'Edit',
      icon: 'undo-2',
      alias: ['undo', 'back'],
      shortcut: { key: 'z', mod: true },
    }),
  ];
}
