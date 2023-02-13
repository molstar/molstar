/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateBuilder, StateObjectRef } from '../../mol-state';
import { StructureFromPetworld } from './model';
import { PetworldColorThemeProvider } from './color';
import { ColorNames } from '../../mol-util/color/names';
import { StateTransforms } from '../../mol-plugin-state/transforms';
import { distinctColors } from '../../mol-util/color/distinct';
import { Color } from '../../mol-util/color';
import { SpacefillRepresentationProvider } from '../../mol-repr/structure/representation/spacefill';
import { StructureRepresentation3D } from '../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../mol-plugin/context';
import { PluginStateObject } from '../../mol-plugin-state/objects';

function getSpacefillParams(color: Color) {
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: true,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels: [
                    { minDistance: 1, maxDistance: 1153, overlap: 0, stride: 1, scaleBias: 1 },
                    { minDistance: 1153, maxDistance: 2307, overlap: 100, stride: 10, scaleBias: 3 },
                    { minDistance: 2307, maxDistance: 10000000, overlap: 300, stride: 50, scaleBias: 2.5 },
                ],
                quality: 'lowest', // avoid 'auto', triggers boundary calc
            },
        },
        colorTheme: {
            name: PetworldColorThemeProvider.name,
            params: {
                value: color,
                saturation: 0,
                lightness: 0,
            }
        },
        sizeTheme: {
            name: 'physical',
            params: {
                scale: 1,
            }
        },
    };
}

export async function createPetworldHierarchy(plugin: PluginContext, trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>) {
    const tr = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory)?.obj?.data;
    if (!tr) return;

    plugin.managers.interactivity.setProps({ granularity: 'chain' });
    plugin.canvas3d?.setProps({
        multiSample: { mode: 'off' },
        cameraClipping: { far: false, minNear: 50 },
        renderer: {
            colorMarker: true,
            highlightColor: Color(0xffffff),
            highlightStrength: 0,
            selectStrength: 0,
            dimColor: Color(0xffffff),
            dimStrength: 0,
            markerPriority: 2,
            interiorColorFlag: false,
            interiorDarkening: 0.15,
        },
        marking: {
            enabled: false,
            highlightEdgeColor: Color(0x394e65),
            selectEdgeStrength: 0,
            ghostEdgeStrength: 1,
            innerEdgeFactor: 2.5,
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

    const state = plugin.state.data;

    const group = await state.build()
        .to(trajectory)
        .group(StateTransforms.Misc.CreateGroup, { label: 'root' }, { tags: 'Entity', state: { isCollapsed: true } })
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
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (let i = 0; i < tr.frameCount; i++) {
                build = build
                    .to(group)
                    .apply(StructureFromPetworld, { modelIndex: i }, { tags: 'Entity' })
                    .apply(StructureRepresentation3D, getSpacefillParams(colors[i]));
            }
            await build.commit();
        } catch (e) {
            console.error(e);
            plugin.log.error(e);
        } finally {
            plugin.animationLoop.start();
        }
    }).run();
}
