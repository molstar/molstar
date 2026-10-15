/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { SegmentRepresentationProvider } from '@molstar/graphics/repr/volume/segment';
import { VolumeSegmentColorThemeProvider } from '@molstar/graphics/theme/color/volume-segment';
import { UniformSizeThemeProvider } from '@molstar/graphics/theme/size/uniform';

/** The segment volume representation with the providers of its default color (volume-segment) and size (uniform) themes. */
export const Segment: PluginRegistryEntry = {
  volume: {
    representations: [SegmentRepresentationProvider],
    themes: { color: [VolumeSegmentColorThemeProvider], size: [UniformSizeThemeProvider] },
  },
};
