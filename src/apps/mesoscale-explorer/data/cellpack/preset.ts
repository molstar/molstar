/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../../../mol-plugin/context';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { StateObjectRef, StateObjectSelector, StateBuilder } from '../../../../mol-state';
import { Color } from '../../../../mol-util/color';
import { distinctColors } from '../../../../mol-util/color/distinct';
import { ColorNames } from '../../../../mol-util/color/names';
import { Hcl } from '../../../../mol-util/color/spaces/hcl';
import { CellpackUniformColorThemeProvider } from './color';
import { CellpackAssembly, EntityStructure } from './model';

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
                    { minDistance: 1, maxDistance: 1153, overlap: 0, stride: 1, scaleBias: 1 },
                    { minDistance: 1153, maxDistance: 2307, overlap: 100, stride: 10, scaleBias: 3 },
                    { minDistance: 2307, maxDistance: 10000000, overlap: 300, stride: 50, scaleBias: 2.5 },
                ],
                quality: 'lowest', // avoid 'auto', triggers boundary calc
            },
        },
        colorTheme: {
            name: CellpackUniformColorThemeProvider.name,
            params: {
                value: color,
                saturation: 0,
                lightness: 0,
            }
        },
        sizeTheme: {
            name: 'physical',
            params: {
                value: 1,
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
