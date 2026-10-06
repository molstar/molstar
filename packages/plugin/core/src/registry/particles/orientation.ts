/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { OrientationParticlesRepresentationProvider } from '@molstar/graphics/repr/particles/representation/orientation';
import { ParticleIndexColorThemeProvider } from '@molstar/graphics/theme/color/particle-index';
import { ParticleSizeThemeProvider } from '@molstar/graphics/theme/size/particle-size';

/** The orientation particles representation with the providers of its default color (particle-index) and size (particle-size) themes. */
export const ParticleOrientation: PluginRegistryEntry = {
  particles: {
    representations: [OrientationParticlesRepresentationProvider],
    themes: { color: [ParticleIndexColorThemeProvider], size: [ParticleSizeThemeProvider] },
  },
};
