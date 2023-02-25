/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateBuilder, StateObjectRef } from '../../../../mol-state';
import { StructureFromPetworld } from './model';
import { PetworldColorThemeProvider } from './color';
import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { distinctColors } from '../../../../mol-util/color/distinct';
import { Color } from '../../../../mol-util/color';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../../../mol-plugin/context';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';

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
                    { minDistance: 1, maxDistance: 1000, overlap: 0, stride: 1, scaleBias: 1 },
                    { minDistance: 1000, maxDistance: 4000, overlap: 500, stride: 10, scaleBias: 3 },
                    { minDistance: 4000, maxDistance: 10000000, overlap: 500, stride: 50, scaleBias: 2.5 },
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
