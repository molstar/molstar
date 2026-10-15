/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { ExternalStructureColorThemeProvider } from './external-structure.js';
import { ExternalVolumeColorThemeProvider } from './external-volume.js';

const ColorThemes = [ExternalStructureColorThemeProvider, ExternalVolumeColorThemeProvider] as const;

/**
 * The `external-structure` and `external-volume` color themes in every scope, as in 5.x. They live in the plugin
 * layer because their parameters list structures and volumes from plugin state.
 */
export const ExternalColorThemes: PluginRegistryEntry = {
  structure: { themes: { color: ColorThemes } },
  volume: { themes: { color: ColorThemes } },
  particles: { themes: { color: ColorThemes } },
};
