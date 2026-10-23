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

/** The trajectory commands (api.md §9, §2 `trajectory`): one per operation, the read aside. */
export function trajectoryCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  const category = 'Models';
  return [
    operationCommand(api, api.trajectory.frame, {
      category,
      icon: 'trajectory',
      title: 'Show a model…',
      alias: ['frame', 'model', 'next model', 'previous model', 'nmr model', 'switch model', 'ensemble', 'trajectory'],
      action: 'Show',
    }),
    operationCommand(api, api.trajectory.play, {
      category,
      icon: 'play',
      title: 'Play the models',
      alias: ['play trajectory', 'animate models', 'ensemble', 'frames'],
    }),
    operationCommand(api, api.trajectory.periodic, {
      category,
      icon: 'trajectory',
      title: 'Show the periodic box…',
      alias: [
        'periodic box',
        'periodic boundary',
        'pbc',
        'make molecules whole',
        'unwrap',
        'wrap into box',
        'center protein in box',
        'no jump',
        'simulation box',
      ],
      action: 'Show',
    }),
    operationCommand(api, api.trajectory.stop, {
      category,
      icon: 'pause',
      title: 'Stop playing the models',
      alias: ['pause trajectory', 'stop trajectory'],
    }),
  ];
}
