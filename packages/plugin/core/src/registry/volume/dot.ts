/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { DotRepresentationProvider } from '@molstar/graphics/repr/volume/dot';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The dot volume representation with the providers of its default color (uniform) and size (uniform) themes. */
export const Dot: PluginRegistryEntry = {
  volume: {
    representations: [DotRepresentationProvider],
    themes: { color: [UniformColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
