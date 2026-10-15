/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';

/** The gaussian-surface structure representation with the providers of its default color (chain-id) and size (physical) themes. */
export const GaussianSurface: PluginRegistryEntry = {
  structure: {
    representations: [GaussianSurfaceRepresentationProvider],
    themes: { color: [ChainIdColorThemeProvider], size: [PhysicalSizeThemeProvider] },
  },
};
