/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { GaussianVolumeRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-volume';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';

/** The gaussian-volume structure representation with the providers of its default color (chain-id) and size (physical) themes. */
export const GaussianVolume: PluginRegistryEntry = {
  structure: {
    representations: [GaussianVolumeRepresentationProvider],
    themes: { color: [ChainIdColorThemeProvider], size: [PhysicalSizeThemeProvider] },
  },
};
