/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { DirectVolumeRepresentationProvider } from '@molstar/graphics/repr/volume/direct-volume';
import { VolumeValueColorThemeProvider } from '@molstar/graphics/theme/color/volume-value';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The direct-volume volume representation with the providers of its default color (volume-value) and size (uniform) themes. */
export const DirectVolume: PluginRegistryEntry = {
  volume: {
    representations: [DirectVolumeRepresentationProvider],
    themes: { color: [VolumeValueColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
