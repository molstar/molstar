/**
 * Copyright (c) 2022-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
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
import { GraphicsMode, MesoscaleGroup, MesoscaleState, getDistinctBaseColors, getGraphicsModeProps, getMesoscaleGroupParams } from '../state';
import { ColorNames } from '../../../../mol-util/color/names';
import { MmcifFormat } from '../../../../mol-model-formats/structure/mmcif';
import { Task } from '../../../../mol-task';

function getSpacefillParams(color: Color, graphics: GraphicsMode) {
    const gmp = getGraphicsModeProps(graphics === 'custom' ? 'quality' : graphics);
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: true,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels: gmp.lodLevels,
                quality: 'lowest', // avoid 'auto', triggers boundary calc
                clip: {
                    variant: 'instance',
                    objects: [],
                },
                clipPrimitive: true,
                approximate: gmp.approximate,
                alphaThickness: gmp.alphaThickness,
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
        // cannot use m.properties.structAsymMap because petworld models
        // may assign the same asymId to multiple entities
        const { label_asym_id, label_entity_id, _rowCount } = m.atomicHierarchy.chains;
        const membraneIds: string[] = [];
        const otherIds: string[] = [];
        const seen = new Set<string>();
        for (let i = 0; i < _rowCount; i ++) {
            const entityId = label_entity_id.value(i);
            if (seen.has(entityId)) continue;

            const asymId = label_asym_id.value(i);
            if (asymId.startsWith('MEM')) {
                membraneIds.push(entityId);
            } else {
                otherIds.push(entityId);
            }
            seen.add(entityId);
        }
        if (membraneIds.length) {
            membrane.push({ modelIndex: i, entityIds: membraneIds });
        }
        if (otherIds.length) {
            other.push({ modelIndex: i, entityIds: otherIds });
        }
    }

    const state = plugin.state.data;
    const graphicsMode = MesoscaleState.get(plugin).graphics;
    const groupParams = getMesoscaleGroupParams(graphicsMode);

    const group = await state.build()
        .toRoot()
        .apply(MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `ent:`, label: 'entity', color: { type: 'generate', illustrative: false, value: ColorNames.white, variability: 20, shift: 0, lightness: 0, alpha: 1, emissive: 0 } }, { tags: ['group:ent:'], state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit({ revertOnError: true });

    await state.build()
        .to(group)
        .apply(MesoscaleGroup, { ...groupParams, index: undefined, tag: `ent:mem`, label: 'Membrane', color: { type: 'uniform', illustrative: false, value: ColorNames.lightgrey, variability: 20, shift: 0, lightness: 0, alpha: 1, emissive: 0 } }, { tags: ['group:ent:mem', 'ent:', '__no_group_color__'], state: { isCollapsed: true, isHidden: groupParams.hidden } })
        .commit();

    const colors = getDistinctBaseColors(other.length, 0);

    await state.transaction(async () => {
        try {
            plugin.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (let i = 0, il = membrane.length; i < il; ++i) {
                build = build
                    .to(cell)
                    .apply(StructureFromPetworld, membrane[i])
                    .apply(StructureRepresentation3D, getSpacefillParams(ColorNames.lightgrey, graphicsMode), { tags: ['ent:mem', '__no_group_color__'] });
            }
            for (let i = 0, il = other.length; i < il; ++i) {
                build = build
                    .to(cell)
                    .apply(StructureFromPetworld, other[i])
                    .apply(StructureRepresentation3D, getSpacefillParams(colors[i], graphicsMode), { tags: ['ent:'] });
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
