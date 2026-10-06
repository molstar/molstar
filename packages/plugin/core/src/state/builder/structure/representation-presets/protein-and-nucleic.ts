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
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { Cartoon } from '@molstar/plugin/registry/structure/cartoon';
import { GaussianSurface } from '@molstar/plugin/registry/structure/gaussian-surface';
import { nucleic, protein } from '@molstar/plugin/state/queries/structure/type';
import { StructureRepresentationPresetProvider, presetSelectionComponent, BuiltInPresetGroupName } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const ProteinAndNucleicPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-protein-and-nucleic',
  alias: 'protein-and-nucleic',
  display: {
    name: 'Protein & Nucleic',
    group: BuiltInPresetGroupName,
    description: 'Shows proteins as Cartoon and RNA/DNA as Gaussian Surface.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      protein: await presetSelectionComponent(plugin, structureCell, protein, 'protein'),
      nucleic: await presetSelectionComponent(plugin, structureCell, nucleic, 'nucleic'),
    };

    const structure = structureCell.obj!.data;
    const gaussianProps = {
      radiusOffset: structure.isCoarseGrained ? 2 : 0,
      smoothness: structure.isCoarseGrained ? 1.0 : 1.5,
    };

    const { update, builder, typeParams, surfaceTypeParams, symmetryColor, symmetryColorParams } = reprBuilder(
      plugin,
      params,
      structure,
    );

    const representations = {
      protein: builder.buildRepresentation(
        update,
        components.protein,
        { type: CartoonRepresentationProvider, typeParams, color: symmetryColor, colorParams: symmetryColorParams },
        { tag: 'protein' },
      ),
      nucleic: builder.buildRepresentation(
        update,
        components.nucleic,
        {
          type: GaussianSurfaceRepresentationProvider,
          typeParams: { ...surfaceTypeParams, ...gaussianProps },
          color: symmetryColor,
          colorParams: symmetryColorParams,
        },
        { tag: 'nucleic' },
      ),
    };

    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name, params.theme?.focus?.params);

    return { components, representations };
  },
});

/** The ProteinAndNucleicPreset preset with the representations and themes it builds. */
export const ProteinAndNucleicPresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { representation: [ProteinAndNucleicPreset] } } },
  Cartoon,
  GaussianSurface,
);
