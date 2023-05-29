/**
 * Copyright (c) 2019-2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Ludovic Autin <ludovic.autin@gmail.com>
 */

import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../../../mol-plugin/context';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { StateObjectRef, StateObjectSelector, StateBuilder } from '../../../../mol-state';
import { Color } from '../../../../mol-util/color';
import { ColorNames } from '../../../../mol-util/color/names';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { MesoscaleExplorerState } from '../../app';
import { MesoscaleGroup, MesoscaleGroupParams, MesoscaleGroupProps, getDistinctBaseColors, getDistinctGroupColors } from '../state';
import { CellpackUniformColorThemeProvider } from './color';
import { CellpackAssembly, EntityStructure } from './model';

type LodLevels = typeof SpacefillRepresentationProvider.defaultValues['lodLevels']

function getSpacefillParams(color: Color, sizeFactor: number, lodLevels: LodLevels) {
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: false,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels,
                quality: 'lowest', // avoid 'auto', triggers boundary calc
                sizeFactor,
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

function getSizeFactor(name: string): number {
    switch (name) {
        case 'dLDL':
            return 2.5;
        case 'iLDL':
            return 5;
        case 'NP_CA':
        case 'POL_CA':
        case 'FactorH1':
        case 'iIgM_Antibody_5mer':
        // case 'MG_271_272_273_274_192MER': // has a coarse and an atomic part
            return 2;
        default: return 1;
    }
}

export async function createCellpackHierarchy(plugin: PluginContext, trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>) {
    const builder = plugin.builders.structure;
    const state = plugin.state.data;
    const customState = plugin.customState as MesoscaleExplorerState;

    const model = await builder.createModel(trajectory, { modelIndex: 0 });
    const entities = model.data!.entities.data;

    const compGroups = new Map<string, StateObjectSelector>();
    const compIds = new Map<string, { idx: number, members: Map<string, number> }>();
    const compColors = new Map<string, Color[]>();

    const funcGroups = new Map<string, StateObjectSelector>();
    const funcIds = new Map<string, { idx: number, size: number }>();
    const funcColors = new Map<string, Color[]>();

    const _groupParams = PD.getDefaultValues(MesoscaleGroupParams);
    const groupParams: MesoscaleGroupProps = {
        ..._groupParams,
        lod: {
            ..._groupParams.lod,
            lodLevels: customState.lodLevels,
        }
    };

    const base = await state.build()
        .to(model)
        .apply(CellpackAssembly, { id: '' })
        .commit();

    const compRoot = await state.build()
        .toRoot()
        .applyOrUpdateTagged('group:comp:', MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `comp:`, label: 'compartment', color: { type: 'custom', value: ColorNames.white, lightness: 0, alpha: 1 } }, { tags: 'group:comp:', state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit();

    const funcRoot = await state.build()
        .toRoot()
        .applyOrUpdateTagged('group:func:', MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `func:`, label: 'function', color: { type: 'custom', value: ColorNames.white, lightness: 0, alpha: 1 } }, { tags: 'group:func:', state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit();

    if (entities._rowCount > 1) {
        for (let i = 0; i < entities._rowCount; i++) {
            const description = entities.pdbx_description.value(i)[0] || 'unkown compartment';
            const d = description.split('.');
            const n = d.slice(0, -1).join('.');
            const l = d.at(-1)!;
            if (!compIds.has(n)) {
                compIds.set(n, { idx: compIds.size, members: new Map() });
            }
            const cm = compIds.get(n)!;
            cm.members.set(l, cm.members.size);

            const f = entities.details.value(i) || 'unknown function';
            if (!funcIds.has(f)) {
                funcIds.set(f, { idx: funcIds.size, size: 0 });
            }
            funcIds.get(f)!.size += 1;
        }

        //

        const baseCompColors = getDistinctBaseColors(compIds.size);
        const compIdEntries = Array.from(compIds.entries());
        for (let i = 0; i < compIdEntries.length; ++i) {
            const [n, m] = compIdEntries[i];
            const groupColors = getDistinctGroupColors(m.members.size, baseCompColors[i]);
            compColors.set(n, groupColors);
        }

        //

        const baseFuncColors = getDistinctBaseColors(funcIds.size);
        const funcIdEntries = Array.from(funcIds.entries());
        for (let i = 0; i < funcIdEntries.length; ++i) {
            const [n, m] = funcIdEntries[i];
            const groupColors = getDistinctGroupColors(m.size, baseFuncColors[i]);
            funcColors.set(n, groupColors);
        }

        //

        for (let i = 0; i < entities._rowCount; i++) {
            const description = entities.pdbx_description.value(i)[0] || 'unkown compartment';
            const nodes = description.split('.');
            for (let j = 0, jl = nodes.length - 1; j < jl; ++j) {
                const n = nodes.slice(0, j + 1).join('.');
                const p = nodes.slice(0, j).join('.');
                if (!compGroups.has(n)) {
                    const colorIdx = compIds.get(n)?.idx;
                    const color = colorIdx !== undefined ? baseCompColors[colorIdx] : ColorNames.white;
                    const label = nodes[j];
                    const parent = compGroups.get(p) ?? compRoot;
                    parent.cell!.state.isCollapsed = false;
                    const group = await state.build()
                        .to(parent)
                        .applyOrUpdateTagged(`group:comp:${n}`, MesoscaleGroup, { ...groupParams, root: parent === compRoot, index: colorIdx, tag: `comp:${n}`, label, color: { type: 'generate', value: color, lightness: 0, alpha: 1 } }, { tags: `comp:${p}`, state: { isCollapsed: true, isHidden: groupParams.hidden } })
                        .commit({ revertOnError: true });
                    compGroups.set(n, group);
                }
            }

            const f = entities.details.value(i) || 'unknown function';
            if (!funcGroups.has(f)) {
                const colorIdx = funcIds.get(f)?.idx;
                const color = colorIdx !== undefined ? baseFuncColors[colorIdx] : ColorNames.white;
                const group = await state.build()
                    .to(funcRoot)
                    .applyOrUpdateTagged(`group:func:${f}`, MesoscaleGroup, { ...groupParams, index: colorIdx, tag: `func:${f}`, label: f, color: { type: 'custom', value: color, lightness: 0, alpha: 1 } }, { tags: 'func:', state: { isCollapsed: true, isHidden: groupParams.hidden } })
                    .commit({ revertOnError: true });
                funcGroups.set(f, group);
            }
        }

        //

        await state.transaction(async () => {
            try {
                plugin.animationLoop.stop({ noDraw: true });
                let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
                for (let i = 0; i < entities._rowCount; i++) {
                    const description = entities.pdbx_description.value(i)[0] || 'model';
                    const d = description.split('.');
                    const n = d.slice(0, -1).join('.');
                    const l = d.at(-1)!;

                    const f = entities.details.value(i) || 'unknown function';

                    const color = compColors.get(n)![compIds.get(n)!.members.get(l)!];
                    const sizeFactor = getSizeFactor(l);

                    build = build
                        .to(base)
                        .apply(EntityStructure, { entityId: entities.id.value(i) })
                        .apply(StructureRepresentation3D, getSpacefillParams(color, sizeFactor, customState.lodLevels), { tags: [`comp:${n}`, `func:${f}`] });
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
            .apply(EntityStructure, { entityId: entities.id.value(0) })
            .apply(StructureRepresentation3D, getSpacefillParams(ColorNames.lightgray, 1, customState.lodLevels), { tags: [`comp:`, `func:`] })
            .commit();
    }
}