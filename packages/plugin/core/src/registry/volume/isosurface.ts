/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The isosurface volume representation with the providers of its default color (uniform) and size (uniform) themes. */
export const Isosurface: PluginRegistryEntry = {
  volume: {
    representations: [IsosurfaceRepresentationProvider],
    themes: { color: [UniformColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
