/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { ParticleTargetRepresentationProvider } from '@molstar/graphics/repr/particles/representation/target/representation';
import { ParticleEntityColorThemeProvider } from '@molstar/graphics/theme/color/particle-entity';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The target particles representation with the providers of its default color (particle-entity) and size (uniform) themes. */
export const ParticleTarget: PluginRegistryEntry = {
  particles: {
    representations: [ParticleTargetRepresentationProvider],
    themes: { color: [ParticleEntityColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
