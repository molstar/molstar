/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { Structure } from '@molstar/model/model/structure';
import { PluginConfig } from '@molstar/plugin/config';
import { assertUnreachable } from '@molstar/core/util/type-helpers';
import { CoarseSurfacePreset } from './coarse-surface.js';
import { PolymerCartoonPreset } from './polymer-cartoon.js';
import { PolymerAndLigandPreset } from './polymer-and-ligand.js';
import { AtomicDetailPreset } from './atomic-detail.js';
import { StructureRepresentationPresetProvider } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;

export const AutoPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-auto',
  alias: 'auto',
  display: {
    name: 'Automatic',
    description:
      'Show representations based on the size of the structure. Smaller structures are shown with more detail than larger ones, ranging from atomistic display to coarse surfaces.',
  },
  params: () => CommonParams,
  apply(ref, params, plugin) {
    const structure = StateObjectRef.resolveAndCheck(plugin.state.data, ref)?.obj?.data;
    if (!structure) return {};

    const thresholds = plugin.config.get(PluginConfig.Structure.SizeThresholds) || Structure.DefaultSizeThresholds;
    const size = Structure.getSize(structure, thresholds);

    const gapFraction = structure.polymerResidueCount / structure.polymerGapCount;

    switch (size) {
      case Structure.Size.Gigantic:
      case Structure.Size.Huge:
        return CoarseSurfacePreset.apply(ref, params, plugin);
      case Structure.Size.Large:
        return PolymerCartoonPreset.apply(ref, params, plugin);
      case Structure.Size.Medium:
        if (gapFraction > 3) {
          return PolymerAndLigandPreset.apply(ref, params, plugin);
        } // else fall through
      case Structure.Size.Small:
        // `showCarbohydrateSymbol: true` is nice, e.g., for PDB 1aga
        return AtomicDetailPreset.apply(ref, { ...params, showCarbohydrateSymbol: true }, plugin);
      default:
        assertUnreachable(size);
    }
  },
});
