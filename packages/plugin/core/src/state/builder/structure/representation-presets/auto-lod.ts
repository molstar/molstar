/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { Vec3 } from '@molstar/core/math/linear-algebra/3d/vec3';
import { mergeRegistryEntries } from '@molstar/plugin/registry/merge';
import type { PluginRegistryEntry } from '@molstar/plugin/spec';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { BallAndStick } from '@molstar/plugin/registry/structure/ball-and-stick';
import { Cartoon } from '@molstar/plugin/registry/structure/cartoon';
import { GaussianSurface } from '@molstar/plugin/registry/structure/gaussian-surface';
import { StructureRepresentationPresetProvider, presetStaticComponent } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const AutoLodPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-auto-lod',
  alias: 'auto-lod',
  display: {
    name: 'Automatic Detail',
    group: 'Miscellaneous',
    description: 'Shows more (or less) detailed representations automatically based on camera distance.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      all: await presetStaticComponent(plugin, structureCell, 'all'),
    };

    const structure = structureCell.obj!.data;

    const {
      update,
      builder,
      typeParams,
      surfaceTypeParams,
      color,
      symmetryColor,
      symmetryColorParams,
      ballAndStickColor,
    } = reprBuilder(plugin, params, structure);

    const representations = {
      gaussianSurface: builder.buildRepresentation(
        update,
        components.all,
        {
          type: GaussianSurfaceRepresentationProvider,
          typeParams: { ...surfaceTypeParams, lod: Vec3.create(30, 10000000, 100) },
          color: symmetryColor,
          colorParams: symmetryColorParams,
        },
        { tag: 'gaussian-surface' },
      ),
      cartoon: builder.buildRepresentation(
        update,
        components.all,
        {
          type: CartoonRepresentationProvider,
          typeParams: { ...typeParams, lod: Vec3.create(-20, 300, 100) },
          color: symmetryColor,
          colorParams: symmetryColorParams,
        },
        { tag: 'cartoon' },
      ),
      ballAndStick: builder.buildRepresentation(
        update,
        components.all,
        {
          type: BallAndStickRepresentationProvider,
          typeParams: { ...typeParams, lod: Vec3.create(-20, 40, 20) },
          color,
          colorParams: ballAndStickColor,
        },
        { tag: 'ball-and-stick' },
      ),
    };

    await update.commit({ revertOnError: false });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name, params.theme?.focus?.params);

    return { components, representations };
  },
});

/** The AutoLodPreset preset with the representations and themes it builds. */
export const AutoLodPresetEntry: PluginRegistryEntry = mergeRegistryEntries(
  { structure: { presets: { representation: [AutoLodPreset] } } },
  GaussianSurface,
  Cartoon,
  BallAndStick,
);
