/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { OrientationRepresentationProvider } from '@molstar/graphics/repr/structure/representation/orientation';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The orientation structure representation with the providers of its default color (chain-id) and size (uniform) themes. */
export const Orientation: PluginRegistryEntry = {
  structure: {
    representations: [OrientationRepresentationProvider],
    themes: { color: [ChainIdColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
