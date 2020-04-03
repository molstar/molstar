/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateObjectRef } from '../../../../mol-state';
import { StructureRepresentationPresetProvider, presetStaticComponent } from '../../../../mol-plugin-state/builder/structure/representation-preset';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';

export const CellpackPackingsPresetParams = {
    traceOnly: PD.Boolean(true),
    polymerOnly: PD.Boolean(true),
    representation: PD.Select('gaussian-surface', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'])),
    hue: PD.Interval([0, 360])
}
export type CellpackPackingsPresetParams = PD.ValuesFor<typeof CellpackPackingsPresetParams>

export const CellpackPackingsPreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-cellpack',
    display: { name: 'CellPack' },
    params: () => CellpackPackingsPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const reprProps = {
            ignoreHydrogens: true,
            traceOnly: params.traceOnly
        };
        const components = {
            structure: params.polymerOnly
                ? await presetStaticComponent(plugin, structureCell, 'polymer')
                : await presetStaticComponent(plugin, structureCell, 'all')
        };
        const selectionType = params.polymerOnly ? 'polymer' : 'all'

        if (params.representation === 'gaussian-surface') {
            Object.assign(reprProps, {
                quality: 'custom', resolution: 10, radiusOffset: 2, doubleSided: false
            })
        } else if (params.representation === 'spacefill' && params.traceOnly) {
            Object.assign(reprProps, { sizeFactor: 2 })
        }

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
        const color = 'model-index'
        const colorParams = {
            palette: {
                name: 'generate',
                params: {
                    hue: params.hue, chroma: [30, 80], luminance: [15, 85],
                    clusteringStepCount: 50, minSampleCount: 800,
                    maxCount: 75
                }
            }
        }
        const representations = {
            structure: builder.buildRepresentation<any>(update, components.structure, { type: params.representation, typeParams: { ...typeParams, ...reprProps }, color, colorParams }, { tag: selectionType })
        };

        await plugin.updateDataState(update, { revertOnError: true });
        return { components, representations };
    }
});