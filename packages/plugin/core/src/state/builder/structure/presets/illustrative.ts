/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { StructureRepresentationPresetProvider, presetStaticComponent } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const IllustrativePreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-illustrative',
  alias: 'illustrative',
  display: {
    name: 'Illustrative',
    group: 'Miscellaneous',
    description: 'Show everything in spacefill representation with illustrative colors and ignore light.',
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

    const { update, builder, typeParams, color } = reprBuilder(plugin, params, structure);

    const representations = {
      all: builder.buildRepresentation(
        update,
        components.all,
        {
          type: 'spacefill',
          typeParams: { ...typeParams, ignoreLight: true },
          color: 'illustrative',
          colorParams: { style: { name: 'entity-id', params: { overrideWater: true } } },
        },
        { tag: 'all' },
      ),
    };
    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name ?? color, params.theme?.focus?.params);

    return { components, representations };
  },
});
