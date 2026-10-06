/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { PuttyRepresentationProvider } from '@molstar/graphics/repr/structure/representation/putty';
import { ChainIdColorThemeProvider } from '@molstar/graphics/theme/color/chain-id';
import { UncertaintySizeThemeProvider } from '@molstar/graphics/theme/size/uncertainty';

/** The putty structure representation with the providers of its default color (chain-id) and size (uncertainty) themes. */
export const Putty: PluginRegistryEntry = {
  structure: {
    representations: [PuttyRepresentationProvider],
    themes: { color: [ChainIdColorThemeProvider], size: [UncertaintySizeThemeProvider] },
  },
};
