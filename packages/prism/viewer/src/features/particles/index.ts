/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParticlesApi, ParticlesPlugin } from './operations.js';
import type { ViewerExtension } from '../../spec.js';
import { particleCommands } from './commands.js';

declare module '../../api.js' {
  // The `particles` area (api.md §2): the particle lists, their representations, visibility and removal.
  interface ExtensionApis {
    particles: ParticlesApi;
  }
}

// The `particles` feature (viewer-spec.md §3.1): Mol*'s particles (`ParticlesRegistry`: the formats, transformers,
// representations and themes) and the `particles` area with its commands. Lazy (`LazyParticles`): a particle file, an
// mmCIF holding particles or a session with particles loads it; the Structure panel shows what is loaded.

export const ParticlesFeature: ViewerExtension = {
  id: 'particles',
  title: 'Particles',
  builtIn: true,
  plugin: ParticlesPlugin,
  api: { particles: (c) => new ParticlesApi(c) },
  commands: (api) => particleCommands(api),
};
