/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParticlesApi } from '../../api/particles.js';
import { VolumeApi } from '../../api/volume.js';
import { operationCommand } from '../../commands/operations.js';
import { volumeCommands } from '../../commands/volume.js';
import type { ViewerExtension } from '../../spec.js';

// The `volume` feature (viewer-spec.md §3.1): maps and particles: the `volume` area (levels, the transfer function,
// streaming) and `particles`, and their commands. Their inspectors are the Structure panel's.

export const VolumeFeature: ViewerExtension = {
  id: 'volume',
  title: 'Volumes',
  builtIn: true,
  api: {
    volume: (c) => new VolumeApi(c),
    particles: (c) => new ParticlesApi(c),
  },
  commands: (api) => [
    ...volumeCommands(api),
    operationCommand(api, api.volume.streamingView, {
      category: 'Volumes',
      icon: 'volume',
      title: 'Set the streamed density view…',
      alias: ['streaming view', 'density around focus', 'density detail', 'whole cell density', 'streaming detail'],
      action: 'Set',
      // The channels' iso values and colors are the streaming row's inspector's (and the tool's).
      values: { channels: {} },
    }),
  ],
};
