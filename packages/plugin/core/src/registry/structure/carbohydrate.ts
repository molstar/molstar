/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { CarbohydrateRepresentationProvider } from '@molstar/graphics/repr/structure/representation/carbohydrate';
import { CarbohydrateSymbolColorThemeProvider } from '@molstar/graphics/theme/color/carbohydrate-symbol';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The carbohydrate structure representation with the providers of its default color (carbohydrate-symbol) and size (uniform) themes. */
export const Carbohydrate: PluginRegistryEntry = {
  structure: {
    representations: [CarbohydrateRepresentationProvider],
    themes: { color: [CarbohydrateSymbolColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
