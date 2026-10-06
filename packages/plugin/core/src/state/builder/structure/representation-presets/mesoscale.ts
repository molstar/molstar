/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Spheres } from '@molstar/graphics/geo/geometry/spheres/spheres';
import { StructureRepresentationPresetProvider, presetStaticComponent } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

type MesoscaleGraphicsMode = keyof typeof Spheres.LodLevelsPresets;
const MesoscaleGraphicsOptions = PD.arrayToOptions(Object.keys(Spheres.LodLevelsPresets) as MesoscaleGraphicsMode[]);
function getMesoscaleLodLevels(mode: MesoscaleGraphicsMode) {
  return Spheres.LodLevelsPresets[mode];
}

export const MesoscalePreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-mesoscale',
  alias: 'mesoscale',
  display: {
    name: 'Mesoscale',
    group: 'Miscellaneous',
    description:
      'Show everything in spacefill representation with instance-granularity and level-of-detail tuned for large particle scenes.',
  },
  params: () => ({
    ...CommonParams,
    graphics: PD.Select<MesoscaleGraphicsMode>('quality', MesoscaleGraphicsOptions),
  }),
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      all: await presetStaticComponent(plugin, structureCell, 'all'),
    };

    const structure = structureCell.obj!.data;

    const { update, builder, typeParams, color } = reprBuilder(plugin, params, structure);

    const graphics: MesoscaleGraphicsMode = params.graphics ?? 'quality';
    const lodLevels = getMesoscaleLodLevels(graphics);
    const approximate = graphics !== 'quality' && graphics !== 'ultra';
    const alphaThickness = graphics === 'performance' ? 15 : 12;

    const representations = {
      all: builder.buildRepresentation(
        update,
        components.all,
        {
          type: 'spacefill',
          typeParams: {
            ...typeParams,
            instanceGranularity: true,
            lodLevels,
            approximate,
            alphaThickness,
            clipPrimitive: true,
          },
          color: color || 'entity-id',
        },
        { tag: 'all' },
      ),
    };

    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name ?? color, params.theme?.focus?.params);

    return { components, representations };
  },
});
