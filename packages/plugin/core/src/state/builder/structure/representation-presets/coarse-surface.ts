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
import { StructureRepresentationPresetProvider, presetStaticComponent, BuiltInPresetGroupName } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const CoarseSurfacePreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-coarse-surface',
  alias: 'coarse-surface',
  display: {
    name: 'Coarse Surface',
    group: BuiltInPresetGroupName,
    description: 'Shows polymers and lipids as coarse Gaussian Surface.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      polymer: await presetStaticComponent(plugin, structureCell, 'polymer'),
      lipid: await presetStaticComponent(plugin, structureCell, 'lipid'),
    };

    const structure = structureCell.obj!.data;
    const thresholds = plugin.config.get(PluginConfig.Structure.SizeThresholds) || Structure.DefaultSizeThresholds;
    const size = Structure.getSize(structure, thresholds);
    const gaussianProps = Object.create(null);
    if (size === Structure.Size.Gigantic) {
      Object.assign(gaussianProps, {
        traceOnly: !structure.isCoarseGrained,
        radiusOffset: 2,
        smoothness: 1,
      });
    } else if (size === Structure.Size.Huge) {
      Object.assign(gaussianProps, {
        radiusOffset: structure.isCoarseGrained ? 2 : 0,
        smoothness: 1,
      });
    } else if (structure.isCoarseGrained) {
      Object.assign(gaussianProps, {
        radiusOffset: 2,
        smoothness: 1,
      });
    }

    const { update, builder, surfaceTypeParams, symmetryColor, symmetryColorParams } = reprBuilder(
      plugin,
      params,
      structure,
    );

    const representations = {
      polymer: builder.buildRepresentation(
        update,
        components.polymer,
        {
          type: 'gaussian-surface',
          typeParams: { ...surfaceTypeParams, ...gaussianProps },
          color: symmetryColor,
          colorParams: symmetryColorParams,
        },
        { tag: 'polymer' },
      ),
      lipid: builder.buildRepresentation(
        update,
        components.lipid,
        {
          type: 'gaussian-surface',
          typeParams: { ...surfaceTypeParams, ...gaussianProps },
          color: symmetryColor,
          colorParams: symmetryColorParams,
        },
        { tag: 'lipid' },
      ),
    };

    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name, params.theme?.focus?.params);

    return { components, representations };
  },
});
