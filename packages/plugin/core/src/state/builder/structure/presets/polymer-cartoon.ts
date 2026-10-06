/**
 * Copyright (c) 2019-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Adam Midlik <midlik@gmail.com>
 * @author Gianluca Tomasello <giagitom@gmail.com>
 */

import { StateObjectRef } from '@molstar/core/state';
import { StructureRepresentationPresetProvider, presetStaticComponent, BuiltInPresetGroupName } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const PolymerCartoonPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-polymer-cartoon',
  alias: 'polymer-cartoon',
  display: {
    name: 'Polymer Cartoon',
    group: BuiltInPresetGroupName,
    description: 'Shows polymers as Cartoon.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      polymer: await presetStaticComponent(plugin, structureCell, 'polymer'),
    };

    const structure = structureCell.obj!.data;

    const { update, builder, typeParams, symmetryColor, symmetryColorParams } = reprBuilder(plugin, params, structure);

    const representations = {
      polymer: builder.buildRepresentation(
        update,
        components.polymer,
        { type: 'cartoon', typeParams, color: symmetryColor, colorParams: symmetryColorParams },
        { tag: 'polymer' },
      ),
    };

    await update.commit({ revertOnError: true });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name, params.theme?.focus?.params);

    return { components, representations };
  },
});
