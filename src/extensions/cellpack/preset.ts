/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { StateBuilder, StateObjectRef, StateObjectSelector } from '../../mol-state';
import { StructureRepresentationPresetProvider, presetStaticComponent } from '../../mol-plugin-state/builder/structure/representation-preset';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { CellPackGenerateColorThemeProvider } from './color/generate';
import { CellPackInfoProvider } from './property';
import { CellPackColorThemeProvider } from './color/basic';
import { EntityStructure, CellpackAssembly } from './state';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { Color } from '../../mol-util/color';
import { distinctColors } from '../../mol-util/color/distinct';
import { Hcl } from '../../mol-util/color/spaces/hcl';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { SpacefillRepresentationProvider } from '../../mol-repr/structure/representation/spacefill';
import { PluginContext } from '../../mol-plugin/context';
import { PluginStateObject } from '../../mol-plugin-state/objects';

export const CellpackPackingPresetParams = {
    traceOnly: PD.Boolean(false),
    ignoreLight: PD.Boolean(true),
    representation: PD.Select('spacefill', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'] as const)),
    uniformColor: PD.Optional(PD.Color(Color(0xFFFFFF))),
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
            let color: string = (structureCell.obj && CellPackInfoProvider.get(structureCell.obj.data).value) ? CellPackGenerateColorThemeProvider.name : CellPackColorThemeProvider.name;
            let colorParams = {};
            if (params.uniformColor !== undefined) {
                color = 'cellpack-uniform';
                colorParams = {
                    value: params.uniformColor,
                    saturation: 0,
                    lightness: 0,
                };
            }

            const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
            const representations = {
                all: builder.buildRepresentation<any>(update, components.all, { type: 'spacefill', typeParams: { ...typeParams, ...reprPropsSpacefill }, color, colorParams }, { tag: 'all' }),
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

function getSpacefillParams(color: Color) {
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: false,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels: [
                    { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1 },
                    { minDistance: 1000, maxDistance: 4000, overlap: 500, stride: 10 },
                    { minDistance: 4000, maxDistance: 10000000, overlap: 500, stride: 50 },
                ],
                quality: 'lowest', // avoid 'auto', triggers boundary calc
            },
        },
        colorTheme: {
            name: 'cellpack-uniform',
            params: {
                value: color,
                saturation: 0,
                lightness: 0,
            }
        },
        sizeTheme: {
            name: 'uniform',
            params: {
                value: 1.7,
            }
        },
    };
}

export async function createCellpackHierarchy(plugin: PluginContext, trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>) {
    const builder = plugin.builders.structure;
    const state = plugin.state.data;

    const model = await builder.createModel(trajectory, { modelIndex: 0 });
    const entities = model.data!.entities.data;
    const groups = new Map<string, StateObjectSelector>();
    const ids = new Map<string, Map<string, number>>();
    const colors = new Map<string, Color[]>();

    const base = await state.build()
        .to(model)
        .apply(CellpackAssembly, { id: '' })
        .commit();

    if (entities._rowCount > 1) {
        for (let i = 0; i < entities._rowCount; i++) {
            const description = entities.pdbx_description.value(i)[0] || 'model';
            const nodes = description.split('.');
            for (let j = 0, jl = nodes.length - 1; j < jl; ++j) {
                const n = nodes.slice(0, j + 1).join('.');
                const p = nodes.slice(0, j).join('.');
                if (!groups.has(n)) {
                    const parent = groups.get(p) ?? base;
                    const group = await state.build()
                        .to(parent)
                        .group(StateTransforms.Misc.CreateGroup, { label: nodes[j] }, { tags: 'Entity', state: { isCollapsed: true } })
                        .commit({ revertOnError: true });
                    groups.set(n, group);
                }
            }
        }

        for (let i = 0; i < entities._rowCount; i++) {
            const description = entities.pdbx_description.value(i)[0] || 'model';
            const d = description.split('.');
            const n = d.slice(0, -1).join('.');
            const l = d.at(-1)!;

            if (!ids.has(n)) {
                ids.set(n, new Map());
            }
            const m = ids.get(n)!;
            m.set(l, m.size);
        }

        const baseColors = distinctColors(ids.size, {
            hue: [1, 360],
            chroma: [40, 70],
            luminance: [15, 85],
            clusteringStepCount: 50,
            minSampleCount: 800,
        });

        const idEntries = Array.from(ids.entries());
        for (let i = 0; i < idEntries.length; ++i) {
            const hcl = Hcl.fromColor(Hcl(), baseColors[i]);
            const hue = [Math.max(1, hcl[0] - 35), Math.min(360, hcl[0] + 35)] as [number, number];
            const [n, m] = idEntries[i];
            const groupColors = distinctColors(m.size, {
                hue,
                chroma: [30, 80],
                luminance: [15, 85],
                clusteringStepCount: 50,
                minSampleCount: 800,
            });
            groupColors[Math.floor(groupColors.length / 2)] = baseColors[i];
            colors.set(n, groupColors);
        }

        await state.transaction(async () => {
            try {
                plugin.animationLoop.stop({ noDraw: true });
                let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
                for (let i = 0; i < entities._rowCount; i++) {
                    const description = entities.pdbx_description.value(i)[0] || 'model';
                    const d = description.split('.');
                    const n = d.slice(0, -1).join('.');
                    const l = d.at(-1)!;

                    const parent = groups.get(n) || base;
                    build = build
                        .to(parent)
                        .apply(EntityStructure, { entityId: entities.id.value(i) }, { tags: 'Entity' })
                        .apply(StructureRepresentation3D, getSpacefillParams(colors.get(n)![ids.get(n)!.get(l)!]));
                }
                await build.commit();
            } catch (e) {
                console.error(e);
                plugin.log.error(e);
            } finally {
                plugin.animationLoop.start();
            }
        }).run();
    } else {
        await state.build()
            .to(base)
            .group(StateTransforms.Misc.CreateGroup, { label: model.obj?.data.label || 'Model' }, { tags: 'Entity', state: { isCollapsed: true } })
            .apply(EntityStructure, { entityId: entities.id.value(0) }, { tags: 'Entity' })
            .apply(StructureRepresentation3D, getSpacefillParams(ColorNames.lightgray))
            .commit();
    }
}
