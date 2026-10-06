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

export const PolymerAndLigandPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-polymer-and-ligand',
  alias: 'polymer-and-ligand',
  display: {
    name: 'Polymer & Ligand',
    group: BuiltInPresetGroupName,
    description:
      'Shows polymers as Cartoon, ligands as Ball & Stick, carbohydrates as 3D-SNFG and water molecules semi-transparent.',
  },
  params: () => CommonParams,
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      polymer: await presetStaticComponent(plugin, structureCell, 'polymer'),
      ligand: await presetStaticComponent(plugin, structureCell, 'ligand'),
      nonStandard: await presetStaticComponent(plugin, structureCell, 'non-standard'),
      branched: await presetStaticComponent(plugin, structureCell, 'branched', { label: 'Carbohydrate' }),
      water: await presetStaticComponent(plugin, structureCell, 'water'),
      ion: await presetStaticComponent(plugin, structureCell, 'ion'),
      lipid: await presetStaticComponent(plugin, structureCell, 'lipid'),
      coarse: await presetStaticComponent(plugin, structureCell, 'coarse'),
    };

    const structure = structureCell.obj!.data;

    // TODO make configurable
    const waterType = (components.water?.obj?.data?.elementCount || 0) > 50_000 ? 'line' : 'ball-and-stick';
    const lipidType = (components.lipid?.obj?.data?.elementCount || 0) > 20_000 ? 'line' : 'ball-and-stick';

    const {
      update,
      builder,
      typeParams,
      color,
      symmetryColor,
      symmetryColorParams,
      globalColorParams,
      ballAndStickColor,
    } = reprBuilder(plugin, params, structure);

    const representations = {
      polymer: builder.buildRepresentation(
        update,
        components.polymer,
        { type: 'cartoon', typeParams, color: symmetryColor, colorParams: symmetryColorParams },
        { tag: 'polymer' },
      ),
      ligand: builder.buildRepresentation(
        update,
        components.ligand,
        { type: 'ball-and-stick', typeParams, color, colorParams: ballAndStickColor },
        { tag: 'ligand' },
      ),
      nonStandard: builder.buildRepresentation(
        update,
        components.nonStandard,
        { type: 'ball-and-stick', typeParams, color, colorParams: ballAndStickColor },
        { tag: 'non-standard' },
      ),
      branchedBallAndStick: builder.buildRepresentation(
        update,
        components.branched,
        { type: 'ball-and-stick', typeParams: { ...typeParams, alpha: 0.3 }, color, colorParams: ballAndStickColor },
        { tag: 'branched-ball-and-stick' },
      ),
      branchedSnfg3d: builder.buildRepresentation(
        update,
        components.branched,
        { type: 'carbohydrate', typeParams, color, colorParams: globalColorParams },
        { tag: 'branched-snfg-3d' },
      ),
      water: builder.buildRepresentation(
        update,
        components.water,
        {
          type: waterType,
          typeParams: {
            ...typeParams,
            alpha: 0.6,
            visuals: waterType === 'line' ? ['intra-bond', 'element-point'] : undefined,
          },
          color,
          colorParams: { carbonColor: { name: 'element-symbol', params: {} }, ...globalColorParams },
        },
        { tag: 'water' },
      ),
      ion: builder.buildRepresentation(
        update,
        components.ion,
        {
          type: 'ball-and-stick',
          typeParams,
          color,
          colorParams: { carbonColor: { name: 'element-symbol', params: {} }, ...globalColorParams },
        },
        { tag: 'ion' },
      ),
      lipid: builder.buildRepresentation(
        update,
        components.lipid,
        {
          type: lipidType,
          typeParams: { ...typeParams, alpha: 0.6, visuals: lipidType === 'line' ? ['intra-bond'] : undefined },
          color,
          colorParams: { carbonColor: { name: 'element-symbol', params: {} }, ...globalColorParams },
        },
        { tag: 'lipid' },
      ),
      coarse: builder.buildRepresentation(
        update,
        components.coarse,
        { type: 'spacefill', typeParams, color: color || 'chain-id', colorParams: globalColorParams },
        { tag: 'coarse' },
      ),
    };

    await update.commit({ revertOnError: false });
    await updateFocusRepr(plugin, structure, params.theme?.focus?.name, params.theme?.focus?.params);

    return { components, representations };
  },
});
