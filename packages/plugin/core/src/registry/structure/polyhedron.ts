/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { PolyhedronRepresentationProvider } from '@molstar/graphics/repr/structure/representation/polyhedron';
import { ElementSymbolColorThemeProvider } from '@molstar/graphics/theme/color/element-symbol';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The polyhedron structure representation with the providers of its default color (element-symbol) and size (uniform) themes. */
export const Polyhedron: PluginRegistryEntry = {
  structure: {
    representations: [PolyhedronRepresentationProvider],
    themes: { color: [ElementSymbolColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
