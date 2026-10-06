/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StructureRepresentationPresetProvider } from './types.js';

export const EmptyPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-empty',
  alias: 'empty',
  display: { name: 'Empty', description: 'Removes all existing representations.' },
  async apply(ref, params, plugin) {
    return {};
  },
});
