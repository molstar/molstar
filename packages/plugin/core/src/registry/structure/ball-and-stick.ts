/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { ElementSymbolColorThemeProvider } from '@molstar/graphics/theme/color/element-symbol';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';

/** The ball-and-stick structure representation with the providers of its default color (element-symbol) and size (physical) themes. */
export const BallAndStick: PluginRegistryEntry = {
  structure: {
    representations: [BallAndStickRepresentationProvider],
    themes: { color: [ElementSymbolColorThemeProvider], size: [PhysicalSizeThemeProvider] },
  },
};
