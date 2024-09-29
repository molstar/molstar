/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { MmcifFormat } from '../../../../mol-model-formats/structure/mmcif';
import { PluginStateObject } from '../../../../mol-plugin-state/objects';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { PluginContext } from '../../../../mol-plugin/context';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { StateObjectRef, StateObjectSelector, StateBuilder } from '../../../../mol-state';
import { Clip } from '../../../../mol-util/clip';
import { Color } from '../../../../mol-util/color';
import { ColorNames } from '../../../../mol-util/color/names';
import { GraphicsMode, MesoscaleGroup, MesoscaleState, getDistinctBaseColors, getDistinctGroupColors, getGraphicsModeProps, getMesoscaleGroupParams } from '../state';
import { MmcifAssembly, MmcifStructure } from './model';

function getSpacefillParams(color: Color, scaleFactor: number, graphics: GraphicsMode, clipVariant: Clip.Variant) {
    const gmp = getGraphicsModeProps(graphics === 'custom' ? 'quality' : graphics);
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: false,
                instanceGranularity: false,
                ignoreLight: true,
                lodLevels: gmp.lodLevels.map(l => {
                    return {
                        ...l,
                        stride: Math.max(1, Math.round(l.stride / Math.pow(scaleFactor, l.scaleBias)))
                    };
                }),
                quality: 'lowest', // avoid 'auto', triggers boundary calc
                clip: {
                    variant: clipVariant,
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
                value: 1,
            }
        },
    };
}

export async function createMmcifHierarchy(plugin: PluginContext, trajectory: StateObjectRef<PluginStateObject.Molecule.Trajectory>) {
    const builder = plugin.builders.structure;
    const state = plugin.state.data;

    const model = await builder.createModel(trajectory, { modelIndex: 0 });
    const { data: entities, subtype } = model.data!.entities;

    const sd = model.data?.sourceData;
    if (MmcifFormat.is(sd)) {
        const pdbId = sd.data.db.struct.entry_id.value(0);
        MesoscaleState.set(plugin, {
            description: sd.data.db.struct.title.value(0),
            link: pdbId ? `https://www.rcsb.org/structure/${pdbId}` : ''
        });
    }

    const spheresAvgRadius = new Map<string, number>();
    if (model.data!.coarseHierarchy.isDefined) {
        const spheresCount = new Map<string, number>();
        const spheresEntity_id = model.data!.coarseHierarchy.spheres.entity_id;
        const spheresRadius = model.data!.coarseConformation.spheres.radius;
        for (let i = 0, il = spheresEntity_id.rowCount; i < il; ++i) {
            const entitiId = spheresEntity_id.value(i);
            const radius = spheresRadius[i];
            if (!spheresCount.has(entitiId)) {
                spheresCount.set(entitiId, 1);
                spheresAvgRadius.set(entitiId, radius);
            } else {
                spheresCount.set(entitiId, spheresCount.get(entitiId)! + 1);
                spheresAvgRadius.set(entitiId, spheresAvgRadius.get(entitiId)! + radius);
            }
        }
        spheresAvgRadius.forEach((v, k) => {
            spheresAvgRadius.set(k, v / spheresCount.get(k)!);
        });
    }

    const entGroups = new Map<string, StateObjectSelector>();
    const entIds = new Map<string, { idx: number, members: Map<number, number> }>();
    const entColors = new Map<string, Color[]>();

    const graphicsMode = MesoscaleState.get(plugin).graphics;
    const groupParams = getMesoscaleGroupParams(graphicsMode);

    const base = await state.build()
        .to(model)
        .apply(MmcifAssembly, { id: '' })
        .commit();

    const units = base.data!.units;
    const willBeMerged = units.length > 1 && units.every(u => u.conformation.operator.isIdentity);
    const clipVariant = willBeMerged ? 'pixel' : 'instance';

    const entRoot = await state.build()
        .toRoot()
        .apply(MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: 'ent:', label: 'entity', color: { type: 'custom', illustrative: false, value: ColorNames.white, variability: 20, shift: 0, lightness: 0, alpha: 1, emissive: 0 } }, { tags: 'group:ent:', state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit();

    const getEntityType = (i: number) => {
        if (entities.type.value(i) === 'water') return 'water' as const;
        return subtype.value(i) || 'unknown type';
    };

    for (let i = 0; i < entities._rowCount; i++) {
        const t = getEntityType(i);
        if (!entIds.has(t)) {
            entIds.set(t, { idx: entIds.size, members: new Map() });
        }
        const cm = entIds.get(t)!;
        cm.members.set(i, cm.members.size);
    }

    //

    const baseEntColors = getDistinctBaseColors(entIds.size, 0);
    const entIdEntries = Array.from(entIds.entries());
    for (let i = 0; i < entIdEntries.length; ++i) {
        const [t, m] = entIdEntries[i];
        const groupColors = getDistinctGroupColors(m.members.size, baseEntColors[i], 20, 0);
        entColors.set(t, groupColors);
    }

    for (let i = 0; i < entities._rowCount; i++) {
        const t = getEntityType(i);
        if (!entGroups.has(t)) {
            const colorIdx = entIds.get(t)?.idx;
            const color = colorIdx !== undefined ? baseEntColors[colorIdx] : ColorNames.white;
            const group = await state.build()
                .to(entRoot)
                .applyOrUpdateTagged(`group:ent:${t}`, MesoscaleGroup, { ...groupParams, index: colorIdx, tag: `ent:${t}`, label: t, color: { type: 'generate', illustrative: false, value: color, variability: 20, shift: 0, lightness: 0, alpha: 1, emissive: 0 } }, { tags: `ent:`, state: { isCollapsed: true, isHidden: groupParams.hidden } })
                .commit({ revertOnError: true });
            entGroups.set(t, group);
        }
    }

    //

    await state.transaction(async () => {
        try {
            const dependsOn = [base.ref];
            plugin.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (let i = 0; i < entities._rowCount; i++) {
                const t = getEntityType(i);
                const color = entColors.get(t)![entIds.get(t)!.members.get(i)!];
                const scaleFactor = spheresAvgRadius.get(entities.id.value(i)) || 1;

                build = build
                    .toRoot()
                    .apply(MmcifStructure, { structureRef: base.ref, entityId: entities.id.value(i) }, { dependsOn })
                    .apply(StructureRepresentation3D, getSpacefillParams(color, scaleFactor, graphicsMode, clipVariant), { tags: [`ent:${t}`] });
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
