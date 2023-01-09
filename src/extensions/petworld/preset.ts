/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { presetStaticComponent, StructureRepresentationPresetProvider } from '../../mol-plugin-state/builder/structure/representation-preset';
import { StateObjectRef } from '../../mol-state';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { TrajectoryHierarchyPresetProvider } from '../../mol-plugin-state/builder/structure/hierarchy-preset';
import { StructureFromPetworld } from './model';
import { PetworldColorThemeParams, PetworldColorThemeProvider } from './color';
import { PluginStateObject } from '../../mol-plugin-state/objects';
import { ColorNames } from '../../mol-util/color/names';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { distinctColors } from '../../mol-util/color/distinct';
import { Color } from '../../mol-util/color';

export const PetworldPresetParams = {
    traceOnly: PD.Boolean(false),
    ignoreLight: PD.Boolean(true),
    representation: PD.Select('spacefill', PD.arrayToOptions(['gaussian-surface', 'spacefill', 'point', 'orientation'] as const)),
    uniformColor: PD.Color(Color(0xFFFFFF)),
};
export type PetworldPresetParams = PD.ValuesFor<typeof PetworldPresetParams>

const PetworldStructurePreset = StructureRepresentationPresetProvider({
    id: 'preset-structure-representation-petworld',
    display: {
        name: 'PetWorld', group: 'Miscellaneous',
        description: '...'
    },
    params: () => PetworldPresetParams,
    async apply(ref, params, plugin) {
        const structureCell = StateObjectRef.resolveAndCheck(plugin.state.data, ref);
        if (!structureCell) return {};

        const reprProps = {
            ignoreHydrogens: true,
            traceOnly: params.traceOnly,
            instanceGranularity: true,
            ignoreLight: params.ignoreLight,
            lodLevels: [
                { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1 },
                { minDistance: 1000, maxDistance: 4000, overlap: 500, stride: 10 },
                { minDistance: 4000, maxDistance: 10000000, overlap: 500, stride: 50 },
            ],
        };
        const components = {
            all: await presetStaticComponent(plugin, structureCell, 'all'),
        };

        const color = PetworldColorThemeProvider.name;
        const colorParams: PD.Values<PetworldColorThemeParams> = {
            value: params.uniformColor
        };

        const { update, builder, typeParams } = StructureRepresentationPresetProvider.reprBuilder(plugin, {});
        const representations = {
            all: builder.buildRepresentation<any>(update, components.all, { type: 'spacefill', typeParams: { ...typeParams, ...reprProps }, color, colorParams }, { tag: 'all' }),
        };

        await update.commit({ revertOnError: true });
        return { components, representations };
    }
});

//

export const PetworldPreset = TrajectoryHierarchyPresetProvider({
    id: 'preset-trajectory-petworld',
    display: {
        name: 'Petworld', group: 'Preset',
        description: 'Shows Petworld structure.'
    },
    isApplicable: o => {
        // TODO: check if petworld cif
        return o.data.frameCount > 1;
    },
    async apply(trajectory, params, plugin) {
        const tr = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory)?.obj?.data;
        if (!tr) return { };

        plugin.managers.interactivity.setProps({ granularity: 'chain' });
        plugin.managers.structure.component.setOptions({
            ...plugin.managers.structure.component.state.options,
            visualQuality: 'custom',
            ignoreLight: true,
            hydrogens: 'hide-all'
        });
        plugin.canvas3d?.setProps({
            multiSample: { mode: 'off' },
            cameraClipping: { far: false },
            renderer: {
                colorMarker: false,
                interiorColorFlag: false,
                interiorDarkening: 0.15,
            },
            marking: {
                enabled: true,
                ghostEdgeStrength: 1,
            },
            postprocessing: {
                occlusion: {
                    name: 'on',
                    params: {
                        samples: 32,
                        radius: 8.5,
                        bias: 1.3,
                        blurKernelSize: 15,
                        resolutionScale: 1,
                    }
                },
                shadow: {
                    name: 'on',
                    params: {
                        bias: 0.6,
                        maxDistance: 80,
                        steps: 3,
                        tolerance: 1.0,
                    }
                },
                outline: {
                    name: 'on',
                    params: {
                        scale: 1,
                        threshold: 0.33,
                        color: ColorNames.black,
                        includeTransparent: true,
                    }
                }
            }
        });

        const builder = plugin.builders.structure;
        const state = plugin.state.data;

        const structures: StateObjectRef<PluginStateObject.Molecule.Structure>[] = [];

        const group = await state.build()
            .to(trajectory)
            .group(StateTransforms.Misc.CreateGroup, { label: 'root' }, { tags: 'Entity' })
            .commit({ revertOnError: true });

        const colors = distinctColors(tr.frameCount, {
            hue: [1, 360],
            chroma: [30, 80],
            luminance: [15, 85],
            clusteringStepCount: 50,
            minSampleCount: 800,
        });

        await state.transaction(async () => {
            try {
                plugin.animationLoop.stop({ noDraw: true });

                for (let i = 0; i < tr.frameCount; i++) {
                    const structure = await state.build()
                        .to(group)
                        .apply(StructureFromPetworld, { modelIndex: i }, { tags: 'Entity', state: { isCollapsed: true } })
                        .commit({ revertOnError: true });

                    structures.push(structure);

                    await builder.representation.applyPreset(structure, PetworldStructurePreset, {
                        traceOnly: false,
                        ignoreLight: true,
                        representation: 'spacefill',
                        uniformColor: colors[i],
                    });
                }
            } catch (e) {
                console.error(e);
                plugin.log.error(e);
            } finally {
                plugin.animationLoop.start();
            }
        }).run();

        return { structures };
    }
});
