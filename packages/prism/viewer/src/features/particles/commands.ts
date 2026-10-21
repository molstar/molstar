/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { CommandDefinition } from '@molstar/prism-ui/core';
import type { ViewerApi } from '../../api/area.js';
import type { ViewerCommandContext } from '../../api/managers.js';
import { operationCommand } from '../../commands/operations.js';

/** The particle commands (api.md §9, §2 `particles`): one per operation, the reads and the inspector's `update` aside. */
export function particleCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  const { particles } = api;
  if (!particles) return [];
  return [
    operationCommand(api, particles.representation, {
      category: 'Particles',
      icon: 'layers',
      title: 'Show particles as…',
      alias: ['particle representation', 'show particles as', 'spacefill particles', 'orientation', 'fibers'],
      action: 'Show',
      values: { params: {} },
    }),
    operationCommand(api, particles.visibility, {
      category: 'Particles',
      icon: 'eye',
      alias: ['show particles', 'hide particles', 'toggle particles'],
      action: 'Set',
    }),
    operationCommand(api, particles.remove, {
      category: 'Particles',
      icon: 'trash-2',
      alias: ['delete particles', 'remove particles'],
      action: 'Remove',
    }),
  ];
}
