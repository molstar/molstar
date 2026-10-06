/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BlobSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/blob-surface';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { PhysicalSizeThemeProvider } from '@molstar/graphics/theme/size/physical';

/** The blob-surface structure representation with the providers of its default color (chain-id) and size (physical) themes. */
export const BlobSurface: PluginRegistryEntry = {
  structure: {
    representations: [BlobSurfaceRepresentationProvider],
    themes: { color: [ChainIdColorThemeProvider], size: [PhysicalSizeThemeProvider] },
  },
};
