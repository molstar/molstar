/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BackboneRepresentationProvider } from '@molstar/graphics/repr/structure/representation/backbone';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The backbone structure representation with the providers of its default color (chain-id) and size (uniform) themes. */
export const Backbone: PluginRegistryEntry = {
  structure: {
    representations: [BackboneRepresentationProvider],
    themes: { color: [ChainIdColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
