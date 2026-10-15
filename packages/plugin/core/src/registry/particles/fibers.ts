/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { FibersRepresentationProvider } from '@molstar/graphics/repr/particles/representation/fibers';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The fibers particles representation with the providers of its default color (uniform) and size (uniform) themes. */
export const ParticleFibers: PluginRegistryEntry = {
  particles: {
    representations: [FibersRepresentationProvider],
    themes: { color: [UniformColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
