/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { SliceRepresentationProvider } from '@molstar/graphics/repr/volume/slice';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The slice volume representation with the providers of its default color (uniform) and size (uniform) themes. */
export const Slice: PluginRegistryEntry = {
  volume: {
    representations: [SliceRepresentationProvider],
    themes: { color: [UniformColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
