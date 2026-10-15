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
 * The volume and particle commands (operations.md §3, §2 `volume`, `particles`): one per operation, the reads and the
 * inspector's `update` aside.
 */
export function volumeCommands(api: ViewerApi): CommandDefinition<ViewerCommandContext>[] {
  const category = 'Volumes';
  return [
    operationCommand(api, api.volume.isovalue, {
      category,
      icon: 'volume',
      alias: ['iso value', 'isovalue', 'contour', 'contour level', 'sigma', 'threshold'],
      action: 'Set',
    }),
    operationCommand(api, api.volume.representation, {
      category,
      icon: 'layers',
      title: 'Show a volume as…',
      alias: ['volume representation', 'isosurface', 'direct volume', 'slice', 'map style', 'show map as'],
      action: 'Show',
      values: { params: {} },
    }),
    operationCommand(api, api.volume.visibility, {
      category,
      icon: 'eye',
      alias: ['show volume', 'hide volume', 'show map', 'hide map', 'toggle volume'],
      action: 'Set',
    }),
    operationCommand(api, api.volume.remove, {
      category,
      icon: 'trash-2',
      alias: ['delete volume', 'remove volume', 'remove map', 'delete map'],
      action: 'Remove',
    }),
    operationCommand(api, api.volume.streaming, {
      category,
      icon: 'volume',
      alias: ['density', 'stream density', 'volume streaming', 'electron density', 'em map', 'download map'],
      action: 'Set',
    }),
    operationCommand(api, api.particles.representation, {
      category: 'Particles',
      icon: 'layers',
      title: 'Show particles as…',
      alias: ['particle representation', 'show particles as', 'spacefill particles', 'orientation', 'fibers'],
      action: 'Show',
      values: { params: {} },
    }),
    operationCommand(api, api.particles.visibility, {
      category: 'Particles',
      icon: 'eye',
      alias: ['show particles', 'hide particles', 'toggle particles'],
      action: 'Set',
    }),
    operationCommand(api, api.particles.remove, {
      category: 'Particles',
      icon: 'trash-2',
      alias: ['delete particles', 'remove particles'],
      action: 'Remove',
    }),
  ];
}
