/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { LabelRepresentationProvider } from '@molstar/graphics/repr/structure/representation/label';
import { UniformColorThemeProvider } from '@molstar/graphics/theme/color/uniform';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';

/** The label structure representation with the providers of its default color (uniform) and size (physical) themes. */
export const Label: PluginRegistryEntry = {
  structure: {
    representations: [LabelRepresentationProvider],
    themes: { color: [UniformColorThemeProvider], size: [PhysicalSizeThemeProvider] },
  },
};
