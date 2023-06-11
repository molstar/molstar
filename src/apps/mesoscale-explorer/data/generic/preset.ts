/**
 * Copyright (c) 2023 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';
import { StateBuilder, StateObjectSelector } from '../../../../mol-state';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';
import { PluginContext } from '../../../../mol-plugin/context';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { Color } from '../../../../mol-util/color';
import { utf8Read } from '../../../../mol-io/common/utf8';
import { Quat, Vec3 } from '../../../../mol-math/linear-algebra';
import { MesoscaleExplorerState } from '../../app';
import { MesoscaleGroup, MesoscaleGroupParams, MesoscaleGroupProps, getDistinctBaseColors, getDistinctGroupColors } from '../state';
import { ColorNames } from '../../../../mol-util/color/names';
import { StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { ReadFile } from '../../../../mol-plugin-state/transforms/data';
import { ModelFromTrajectory, TrajectoryFromGRO } from '../../../../mol-plugin-state/transforms/model';
import { Euler } from '../../../../mol-math/linear-algebra/3d/euler';
import { Asset } from '../../../../mol-util/assets';
import { Clip } from '../../../../mol-util/clip';
import { StructureFromGeneric } from './model';

type LodLevels = typeof SpacefillRepresentationProvider.defaultValues['lodLevels']

function getSpacefillParams(color: Color, sizeFactor: number, lodLevels: LodLevels, clipVariant: Clip.Variant) {
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: true,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels: lodLevels.map(l => {
                    return {
                        ...l,
                        stride: Math.max(1, Math.round(l.stride / Math.pow(sizeFactor, l.scaleBias)))
                    };
                }),
                quality: 'lowest', // avoid 'auto', triggers boundary calc
                sizeFactor,
                clip: {
                    variant: clipVariant,
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

export async function createGenericHierarchy(ctx: PluginContext, file: Asset.File) {
    const asset = await ctx.runTask(ctx.managers.asset.resolve(file, 'zip'));
    const d = asset.data['instanced_structure.json'];
    const t = utf8Read(d, 0, d.length);
    const json = JSON.parse(t) as { model: string, positions: Vec3[], rotations: Vec3[], function: string }[];
    // console.log(json)

    const funcGroups = new Map<string, StateObjectSelector>();
    const funcIds = new Map<string, { idx: number, members: Map<string, number> }>();
    const funcColors = new Map<string, Color[]>();

    const modelTransforms: { name: string, transforms: Mat4[], function: string }[] = [];
    for (const v of json) {
        const transforms: Mat4[] = [];
        const positions = v.positions || (v as any).position;
        const rotations = v.rotations || (v as any).rotation;
        for (let i = 0, il = positions.length; i < il; ++i) {
            let p = positions[i];
            if (Array.isArray(p[0])) p = (p[0] as any);
            let r = rotations[i];
            if (Array.isArray(r[0])) r = (r[0] as any);
            const e = Euler.fromVec3(Euler(), r);
            const q = Quat.fromEuler(Quat(), e, 'XYZ');
            const m = Mat4.fromQuat(Mat4(), q);
            Vec3.scale(p, p, 10);
            Mat4.setTranslation(m, p);
            transforms.push(m);
        }
        modelTransforms.push({ name: v.model, transforms, function: v.function });

        if (!funcIds.has(v.function)) {
            funcIds.set(v.function, { idx: funcIds.size, members: new Map() });
        }
        const fm = funcIds.get(v.function)!;
        fm.members.set(v.model, fm.members.size);
    }
    // console.log(modelTransforms)

    const state = ctx.state.data;
    const customState = ctx.customState as MesoscaleExplorerState;

    const _groupParams = PD.getDefaultValues(MesoscaleGroupParams);
    const groupParams: MesoscaleGroupProps = {
        ..._groupParams,
        lod: {
            ..._groupParams.lod,
            lodLevels: customState.lodLevels,
        }
    };

    const funcRoot = await state.build()
        .toRoot()
        .applyOrUpdateTagged('group:func:', MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `func:`, label: 'function', color: { type: 'custom', value: ColorNames.white, variablity: 35, lightness: 0, alpha: 1 } }, { tags: 'group:func:', state: { isCollapsed: false, isHidden: groupParams.hidden } })
        .commit();

    const baseFuncColors = getDistinctBaseColors(funcIds.size);
    const funcIdEntries = Array.from(funcIds.entries());
    for (let i = 0; i < funcIdEntries.length; ++i) {
        const [n, m] = funcIdEntries[i];
        const groupColors = getDistinctGroupColors(m.members.size, baseFuncColors[i], 35);
        funcColors.set(n, groupColors);
    }

    for (const { function: f } of modelTransforms) {
        if (!funcGroups.has(f)) {
            const colorIdx = funcIds.get(f)?.idx;
            const color = colorIdx !== undefined ? baseFuncColors[colorIdx] : ColorNames.white;
            const group = await state.build()
                .to(funcRoot)
                .applyOrUpdateTagged(`group:func:${f}`, MesoscaleGroup, { ...groupParams, index: colorIdx, tag: `func:${f}`, label: f, color: { type: 'custom', value: color, variablity: 35, lightness: 0, alpha: 1 } }, { tags: 'func:', state: { isCollapsed: true, isHidden: groupParams.hidden } })
                .commit({ revertOnError: true });
            funcGroups.set(f, group);
        }
    }

    await state.transaction(async () => {
        try {
            ctx.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (const { name, transforms, function: f } of modelTransforms) {
                const d = asset.data[name];
                const t = utf8Read(d, 0, d.length);
                const file = Asset.File(new File([t], name));

                const color = funcColors.get(f)![funcIds.get(f)!.members.get(name)!];
                const clipVariant = transforms.length === 1 ? 'pixel' : 'instance';

                build = build
                    .toRoot()
                    .apply(ReadFile, { file, label: name })
                    .apply(TrajectoryFromGRO)
                    .apply(ModelFromTrajectory, { modelIndex: 0 })
                    .apply(StructureFromGeneric, { transforms, label: name.split('.')[0] })
                    .apply(StructureRepresentation3D, getSpacefillParams(color, 1.5, customState.lodLevels, clipVariant), { tags: [`func:${f}`] });
            }
            await build.commit();
        } catch (e) {
            console.error(e);
            ctx.log.error(e);
        } finally {
            ctx.animationLoop.start();
        }
    }).run();
}
