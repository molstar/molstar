/**
 * Copyright (c) 2022-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { StateBuilder, StateObjectRef } from '../../../../mol-state';
import { StructureFromPetworld } from './model';
import { PetworldColorThemeProvider } from './color';
import { StateTransforms } from '../../../../mol-plugin-state/transforms';
import { Color } from '../../../../mol-util/color';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../../../mol-plugin/context';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { MesoscaleExplorerState } from '../../app';
import { MesoscaleGroup, MesoscaleGroupParams, MesoscaleGroupProps, getDistinctBaseColors } from '../state';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { ColorNames } from '../../../../mol-util/color/names';

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
        .to(trajectory)
        .group(StateTransforms.Misc.CreateGroup, { label: 'root' }, { tags: 'Entity', state: { isCollapsed: true } })
        .apply(MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `ent:`, label: 'entity', color: { type: 'generate', value: ColorNames.white, lightness: 0, alpha: 1 } }, { tags: '', state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit({ revertOnError: true });

    const colors = getDistinctBaseColors(tr.frameCount);

    await state.transaction(async () => {
        try {
            plugin.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (let i = 0; i < tr.frameCount; i++) {
                build = build
                    .to(group)
                    .apply(StructureFromPetworld, { modelIndex: i }, { tags: 'Entity' })
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
