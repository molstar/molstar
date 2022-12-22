/**
 * Copyright (c) 2019-2022 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { StateObjectRef } from '../../mol-state';
import { StructureRepresentationPresetProvider, presetStaticComponent } from '../../mol-plugin-state/builder/structure/representation-preset';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { CellPackGenerateColorThemeProvider } from './color/generate';
import { TrajectoryHierarchyPresetProvider } from '../../mol-plugin-state/builder/structure/hierarchy-preset';
import { CellPackInfoProvider } from './property';
import { CellPackColorThemeProvider } from './color/basic';
import { EntityStructure } from './state';
import { PluginStateObject } from '../../mol-plugin-state/objects';

export const CellpackPackingPresetParams = {
    traceOnly: PD.Boolean(false),
    ignoreLight: PD.Boolean(true),
    representation: PD.Select('spacefill', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'] as const)),
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
            traceOnly: params.traceOnly,
            instanceGranularity: true,
            ignoreLight: params.ignoreLight,
        };

        if (params.representation === 'spacefill') {
            const components = {
                all: await presetStaticComponent(plugin, structureCell, 'all')
            };

            const reprPropsSpacefill = {
                ...reprProps,
                lodLevels: [
                    { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1 },
                    { minDistance: 1000, maxDistance: 4000, overlap: 500, stride: 10 },
                    { minDistance: 4000, maxDistance: 10000000, overlap: 500, stride: 50 },
                ],
            };

            // default is generated
            const color = (structureCell.obj && CellPackInfoProvider.get(structureCell.obj.data).value) ? CellPackGenerateColorThemeProvider.name : CellPackColorThemeProvider.name;

            const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
            const representations = {
                all: builder.buildRepresentation<any>(update, components.all, { type: 'spacefill', typeParams: { ...typeParams, ...reprPropsSpacefill }, color }, { tag: 'all' }),
            };

            await update.commit({ revertOnError: true });
            return { components, representations };
        } else {
            const components = {
                polymer: await presetStaticComponent(plugin, structureCell, 'polymer')
            };

            if (params.representation === 'gaussian-surface') {
                Object.assign(reprProps, {
                    quality: 'custom', resolution: 10, radiusOffset: 2, doubleSided: false
                });
            }

            // default is generated
            const color = CellPackGenerateColorThemeProvider.name;

            const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
            const representations = {
                polymer: builder.buildRepresentation<any>(update, components.polymer, { type: params.representation, typeParams: { ...typeParams, ...reprProps }, color }, { tag: 'polymer' })
            };

            await update.commit({ revertOnError: true });
            return { components, representations };
        }
    }
});

//

export const CellpackMembranePresetParams = {
    ignoreLight: PD.Boolean(false),
    representation: PD.Select('spacefill', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'] as const)),
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
            instanceGranularity: true,
            ignoreLight: params.ignoreLight,
        };
        const components = {
            membrane: await presetStaticComponent(plugin, structureCell, 'all', { label: 'Membrane' })
        };

        if (params.representation === 'spacefill') {
            const reprPropsSpacefill = {
                ...reprProps,
                lodLevels: [
                    { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1 },
                    { minDistance: 1000, maxDistance: 4000, overlap: 500, stride: 10 },
                    { minDistance: 4000, maxDistance: 10000000, overlap: 500, stride: 50 },
                ],
            };

            const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
            const representations = {
                membrane: builder.buildRepresentation<any>(update, components.membrane, { type: 'spacefill', typeParams: { ...typeParams, ...reprPropsSpacefill }, color: 'uniform', colorParams: { value: ColorNames.lightgrey } }, { tag: 'all' }),
            };

            await update.commit({ revertOnError: true });
            return { components, representations };
        } else {
            if (params.representation === 'gaussian-surface') {
                Object.assign(reprProps, {
                    quality: 'custom', resolution: 10, radiusOffset: 2, doubleSided: false
                });
            }

            const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
            const representations = {
                membrane: builder.buildRepresentation(update, components.membrane, { type: params.representation, typeParams: { ...typeParams, ...reprProps }, color: 'uniform', colorParams: { value: ColorNames.lightgrey } }, { tag: 'all' })
            };

            await update.commit({ revertOnError: true });
            return { components, representations };
        }
    }
});

//

export const CellpackPreset = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-default',
    display: {
        name: 'Cellpack', group: 'Preset',
        description: 'Shows Cellpack structure.'
    },
    isApplicable: o => {
        // TODO: check if cellpack cif
        return true;
    },
    async apply(trajectory, params, plugin) {
        const builder = plugin.builders.structure;
        const state = plugin.state.data;

        const model = await builder.createModel(trajectory, { modelIndex: 0 });
        const modelProperties = await builder.insertModelProperties(model);
        const entities = model.data!.entities.data;

        const base = await builder.createStructure(modelProperties || model, { name: 'assembly', params: { id: '1' } }, undefined, entities._rowCount === 1 ? 'Entity' : undefined);

        const structures: StateObjectRef<PluginStateObject.Molecule.Structure>[] = [];

        if (entities._rowCount > 1) {
            for (let i = 0; i < entities._rowCount; ++i) {
                const structure = await state.build()
                    .to(base)
                    .apply(EntityStructure, { entityId: entities.id.value(i) }, { tags: 'Entity' })
                    .commit({ revertOnError: true });

                structures.push(structure);

                await builder.representation.applyPreset(structure, CellpackPackingPreset);
            }
        } else {
            structures.push(base);
            await builder.representation.applyPreset(base, CellpackPackingPreset);
        }

        return { structures };
    }
});
