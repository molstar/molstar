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
import { IndexPairBonds } from '@molstar/model/formats/structure/property/bonds/index-pair';
import { StructConn } from '@molstar/model/formats/structure/property/bonds/struct_conn';
import type { StructureRepresentationRegistry } from '@molstar/graphics/repr/structure/registry';
import { StructureRepresentationPresetProvider, presetStaticComponent, BuiltInPresetGroupName } from './types.js';

const CommonParams = StructureRepresentationPresetProvider.CommonParams;
const reprBuilder = StructureRepresentationPresetProvider.reprBuilder;
const updateFocusRepr = StructureRepresentationPresetProvider.updateFocusRepr;

export const AtomicDetailPreset = StructureRepresentationPresetProvider({
  id: 'preset-structure-representation-atomic-detail',
  alias: 'atomic-detail',
  display: {
    name: 'Atomic Detail',
    group: BuiltInPresetGroupName,
    description: 'Shows everything in atomic detail.',
  },
  params: () => ({
    ...CommonParams,
    showCarbohydrateSymbol: PD.Boolean(false),
  }),
  async apply(ref, params, plugin) {
    const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
    if (!structureCell) return {};

    const components = {
      all: await presetStaticComponent(plugin, structureCell, 'all'),
      branched: undefined,
    };

    const structure = structureCell.obj!.data;
    const highElementCount = structure.elementCount > 100_000; // TODO make configurable
    const veryHighElementCount = structure.elementCount > 1_000_000; // TODO make configurable
    const highUnitCount = structure.units.length > 5_000; // TODO make configurable
    const lowResidueElementRatio =
      structure.atomicResidueCount &&
      structure.elementCount > 1000 &&
      structure.atomicResidueCount / structure.elementCount < 3;

    const m = structure.models[0];
    const bondsGiven = !!IndexPairBonds.Provider.get(m) || StructConn.isExhaustive(m);

    let atomicType: StructureRepresentationRegistry.BuiltIn = 'ball-and-stick';
    if (structure.isCoarseGrained || highUnitCount) {
      atomicType = veryHighElementCount ? 'point' : 'spacefill';
    } else if (lowResidueElementRatio && !bondsGiven) {
      atomicType = 'spacefill';
    } else if (highElementCount) {
      atomicType = 'line';
    }
    const showCarbohydrateSymbol = params.showCarbohydrateSymbol && !highElementCount && !lowResidueElementRatio;

    if (showCarbohydrateSymbol) {
      Object.assign(components, {
        branched: await presetStaticComponent(plugin, structureCell, 'branched', { label: 'Carbohydrate' }),
      });
    }

    const { update, builder, typeParams, color, ballAndStickColor, globalColorParams } = reprBuilder(
      plugin,
      params,
      structure,
    );
    const colorParams =
      lowResidueElementRatio && !bondsGiven
        ? { carbonColor: { name: 'element-symbol', params: {} }, ...globalColorParams }
        : ballAndStickColor;

    const representations = {
      all: builder.buildRepresentation(
        update,
        components.all,
        { type: atomicType, typeParams, color, colorParams },
        { tag: 'all' },
      ),
    };
    if (showCarbohydrateSymbol) {
      Object.assign(representations, {
        snfg3d: builder.buildRepresentation(
          update,
          components.branched,
          {
            type: 'carbohydrate',
            typeParams: { ...typeParams, alpha: 0.4, visuals: ['carbohydrate-symbol'] },
            color,
            colorParams: globalColorParams,
          },
          { tag: 'snfg-3d' },
        ),
      });
    }

    await update.commit({ revertOnError: true });
    await updateFocusRepr(
      plugin,
      structure,
      params.theme?.focus?.name ?? color,
      params.theme?.focus?.params ?? colorParams,
    );

    return { components, representations };
  },
});
