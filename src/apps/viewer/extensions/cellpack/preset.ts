/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateObjectRef } from '../../../../mol-state';
import { StructureRepresentationPresetProvider, presetStaticComponent, presetSelectionComponent } from '../../../../mol-plugin-state/builder/structure/representation-preset';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';

export const CellpackPackingsPresetParams = {
    traceOnly: PD.Boolean(true),
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

        const reprProps = Object.create(null);
        const components = Object.create(null);

        let selectionType = 'polymer'

        if (params.traceOnly) {
            selectionType = 'trace'
            components.polymer = await presetSelectionComponent(plugin, structureCell, 'trace')
        } else {
            components.polymer = await presetStaticComponent(plugin, structureCell, 'polymer')
        }

        if (params.representation === 'gaussian-surface') {
            Object.assign(reprProps, {
                quality: 'custom', resolution: 10, radiusOffset: 2,
                alpha: 1.0, flatShaded: false, doubleSided: false,
                ignoreHydrogens: true
            })
        } else if (params.representation === 'spacefill') {
            if (params.traceOnly) {
                Object.assign(reprProps, { sizeFactor: 2, ignoreHydrogens: true })
            } else {
                Object.assign(reprProps, { ignoreHydrogens: true })
            }
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
            polymer: builder.buildRepresentation<any>(update, components.polymer, { type: params.representation, typeParams: { ...typeParams, ...reprProps }, color, colorParams }, { tag: selectionType })
        };

        await plugin.updateDataState(update, { revertOnError: true });
        return { components, representations };
    }
});