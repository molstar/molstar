/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateBuilder, StateObjectRef } from '../../../../mol-state';
import { StructureFromPetworld } from './model';
import { Color } from '../../../../mol-util/color';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../../../mol-plugin/context';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { MesoscaleExplorerState } from '../../app';
import { MesoscaleGroup, MesoscaleGroupParams, MesoscaleGroupProps, getDistinctBaseColors } from '../state';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { ColorNames } from '../../../../mol-util/color/names';
import { MmcifFormat } from '../../../../mol-model-formats/structure/mmcif';
import { Task } from '../../../../mol-task';

type LodLevels = typeof SpacefillRepresentationProvider.defaultValues['lodLevels']

function getSpacefillParams(color: Color, lodLevels: LodLevels) {
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: true,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels,
                quality: 'lowest', // avoid 'auto', triggers boundary calc
                clip: {
                    variant: 'instance',
                    objects: [],
                },
            },
        },
        colorTheme: {
            name: 'uniform',
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
    const cell = StateObjectRef.resolveAndCheck(plugin.state.data, trajectory);
    const tr = cell?.obj?.data;
    if (!cell || !tr) return;

    if (!MmcifFormat.is(tr.representative.sourceData)) return;

    const membrane: { modelIndex: number, entityIds: string[] }[] = [];
    const other: { modelIndex: number, entityIds: string[] }[] = [];
    for (let i = 0; i < tr.frameCount; ++i) {
        const m = await Task.resolveInContext(tr.getFrameAtIndex(i));
        const membraneIds: string[] = [];
        const otherIds: string[] = [];
        m.properties.structAsymMap.forEach((v, k) => {
            if (k.startsWith('MEM')) {
                membraneIds.push(v.entity_id);
            } else {
                otherIds.push(v.entity_id);
            }
        });
        if (membraneIds.length) {
            membrane.push({ modelIndex: i, entityIds: membraneIds });
        }
        if (otherIds.length) {
            other.push({ modelIndex: i, entityIds: otherIds });
        }
    }

    const state = plugin.state.data;
    const customState = plugin.customState as MesoscaleExplorerState;

    const _groupParams = PD.getDefaultValues(MesoscaleGroupParams);
    const groupParams: MesoscaleGroupProps = {
        ..._groupParams,
        lod: {
            ..._groupParams.lod,
            lodLevels: customState.lodLevels,
        }
    };

    const group = await state.build()
        .toRoot()
        .applyOrUpdateTagged('group:ent:', MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `ent:`, label: 'entity', color: { type: 'generate', value: ColorNames.white, variablity: 35, lightness: 0, alpha: 1 } }, { tags: 'group:ent:', state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit({ revertOnError: true });

    await state.build()
        .to(group)
        .applyOrUpdateTagged(`group:ent:mem`, MesoscaleGroup, { ...groupParams, index: undefined, tag: `ent:mem`, label: 'Membrane', color: { type: 'uniform', value: ColorNames.lightgrey, variablity: 35, lightness: 0, alpha: 1 } }, { tags: `ent:`, state: { isCollapsed: true, isHidden: groupParams.hidden } })
        .commit();

    const colors = getDistinctBaseColors(other.length);

    await state.transaction(async () => {
        try {
            plugin.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (let i = 0, il = membrane.length; i < il; ++i) {
                build = build
                    .to(cell)
                    .apply(StructureFromPetworld, membrane[i])
                    .apply(StructureRepresentation3D, getSpacefillParams(ColorNames.lightgrey, customState.lodLevels), { tags: [`ent:mem`] });
            }
            for (let i = 0, il = other.length; i < il; ++i) {
                build = build
                    .to(cell)
                    .apply(StructureFromPetworld, other[i])
                    .apply(StructureRepresentation3D, getSpacefillParams(colors[i], customState.lodLevels), { tags: [`ent:`] });
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
