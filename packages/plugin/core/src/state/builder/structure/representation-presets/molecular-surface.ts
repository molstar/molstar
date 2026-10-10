/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { MolecularSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/molecular-surface';
import { EntityIdColorThemeProvider } from '@molstar/graphics/theme/color/entity-id';
import { MolecularSurface } from '@molstar/plugin/registry/structure/molecular-surface';
import { StructureRepresentationPresetProvider, presetStaticComponent } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const MolecularSurfacePreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-molecular-surface',
  alias: 'molecular-surface',
  display: {
    name: 'Molecular Surface',
    group: 'Miscellaneous',
    description: 'Show everything in molecular surface representation with illustrative colors.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      all: await presetStaticComponent(plugin, structureCell, 'all'),
      branched: undefined,
    };

    const structure = structureCell.obj!.data;

    const { update, builder, surfaceTypeParams, color } = reprBuilder(plugin, params, structure);

    const representations = {
      all: builder.buildRepresentation(
        update,
        components.all,
        {
          type: MolecularSurfaceRepresentationProvider,
          typeParams: surfaceTypeParams,
          color: EntityIdColorThemeProvider,
          colorParams: { overrideWater: true },
        },
        { tag: 'all' },
      ),
    };
    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name ?? color, params.theme?.focus?.params);

    return { components, representations };
  },
});

/** The MolecularSurfacePreset preset with the representations and themes it builds. */
export const MolecularSurfacePresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { representation: [MolecularSurfacePreset] } } },
  MolecularSurface,
  { structure: { themes: { color: [EntityIdColorThemeProvider] } } },
);
