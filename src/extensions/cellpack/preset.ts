/**
 * Copyright (c) 2020 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateObjectRef } from '../../mol-state';
import { StructureRepresentationPresetProvider, presetStaticComponent } from '../../mol-plugin-state/builder/structure/representation-preset';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { CellPackGenerateColorThemeProvider } from './color/generate';
import { CellPackInfoProvider } from './property';
import { CellPackProvidedColorThemeProvider } from './color/provided';

export const CellpackPackingPresetParams = {
    traceOnly: PD.Boolean(true),
    representation: PD.Select('gaussian-surface', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'])),
};
export type CellpackPackingPresetParams = PD.ValuesFor<typeof CellpackPackingPresetParams>

export const CellpackPackingPreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-cellpack-packing',
    display: { name: 'CellPack Packing' },
    params: () => CellpackPackingPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const reprProps = {
            ignoreHydrogens: true,
            traceOnly: params.traceOnly
        };
        const components = {
            polymer: await presetStaticComponent(plugin, structureCell, 'polymer')
        };

        if (params.representation === 'gaussian-surface') {
            Object.assign(reprProps, {
                quality: 'custom', resolution: 10, radiusOffset: 2, doubleSided: false
            });
        } else if (params.representation === 'spacefill' && params.traceOnly) {
            Object.assign(reprProps, { sizeFactor: 2 });
        }

        const info = structureCell.obj?.data && CellPackInfoProvider.get(structureCell.obj?.data).value;
        const color = info?.colors ? CellPackProvidedColorThemeProvider.name : CellPackGenerateColorThemeProvider.name;

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
        const representations = {
            polymer: builder.buildRepresentation<any>(update, components.polymer, { type: params.representation, typeParams: { ...typeParams, ...reprProps }, color }, { tag: 'polymer' })
        };

        await update.commit({ revertOnError: true });
        return { components, representations };
    }
});

//

export const CellpackMembranePresetParams = {
    representation: PD.Select('gaussian-surface', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'])),
};
export type CellpackMembranePresetParams = PD.ValuesFor<typeof CellpackMembranePresetParams>

export const CellpackMembranePreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-cellpack-membrane',
    display: { name: 'CellPack Membrane' },
    params: () => CellpackMembranePresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const reprProps = {
            ignoreHydrogens: true,
        };
        const components = {
            membrane: await presetStaticComponent(plugin, structureCell, 'all', { label: 'Membrane' })
        };

        if (params.representation === 'gaussian-surface') {
            Object.assign(reprProps, {
                quality: 'custom', resolution: 10, radiusOffset: 2, doubleSided: false
            });
        }

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
        const representations = {
            membrane: builder.buildRepresentation(update, components.membrane, { type: 'gaussian-surface', typeParams: { ...typeParams, ...reprProps }, color: 'uniform', colorParams: { value: ColorNames.lightgrey } }, { tag: 'all' })
        };

        await update.commit({ revertOnError: true });

        return { components, representations };
    }
});