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
import { MesoscaleGroup, MesoscaleGroupParams, MesoscaleGroupProps, updateColors } from '../state';
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
    let manifest: GenericManifest;
    // TODO: remove special handling for martini prototype
    if (asset.data['instanced_structure.json']) {
        const d = asset.data['instanced_structure.json'];
        const t = utf8Read(d, 0, d.length);
        const martini = JSON.parse(t) as { model: string, positions: Vec3[], rotations: Vec3[], function: string }[];
        console.log(martini);
        manifest = martiniToGeneric(martini);
    } else if (asset.data['manifest.json']) {
        const d = asset.data['manifest.json'];
        const t = utf8Read(d, 0, d.length);
        manifest = JSON.parse(t) as GenericManifest;
    } else {
        throw new Error('no manifest found');
    }
    console.log(manifest);

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

    async function addGroup(g: GenericGroup, cell: StateObjectSelector, parent: string) {
        const group = await state.build()
            .to(cell)
            .applyOrUpdateTagged(`group:${g.root}:${g.id}`, MesoscaleGroup, { ...groupParams, index: undefined, tag: `${g.root}:${g.id}`, label: g.label || g.id, color: { type: 'custom', value: ColorNames.white, variablity: 35, lightness: 0, alpha: 1 } }, { tags: g.root === parent ? `${g.root}:` : `${g.root}:${parent}`, state: { isCollapsed: true, isHidden: groupParams.hidden } })
            .commit();

        if (g.children) {
            for (const c of g.children) {
                await addGroup(c, group, g.id);
            }
        }
    }

    for (const r of manifest.roots) {
        const root = await state.build()
            .toRoot()
            .applyOrUpdateTagged(`group:${r.id}:`, MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `${r.id}:`, label: r.label || r.id, color: { type: 'custom', value: ColorNames.white, variablity: 35, lightness: 0, alpha: 1 } }, { tags: `group:${r.id}:`, state: { isCollapsed: false, isHidden: groupParams.hidden } })
            .commit();

        if (r.children) {
            for (const c of r.children!) {
                await addGroup(c, root, r.id);
            }
        }
    }

    const p = Vec3();
    const r = Vec3();

    await state.transaction(async () => {
        try {
            ctx.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (const e of manifest.entities) {
                const d = asset.data[e.file];
                const t = utf8Read(d, 0, d.length);
                const file = Asset.File(new File([t], e.file));

                const { positions, rotations } = e.instances;
                const transforms: Mat4[] = [];
                for (let i = 0, il = positions.length; i < il; i += 3) {
                    Vec3.fromArray(p, positions, i);
                    Vec3.fromArray(r, rotations, i);
                    const e = Euler.fromVec3(Euler(), r);
                    const q = Quat.fromEuler(Quat(), e, 'XYZ');
                    const m = Mat4.fromQuat(Mat4(), q);
                    Mat4.setTranslation(m, p);
                    transforms.push(m);
                }

                const color = ColorNames.skyblue;
                const clipVariant = transforms.length === 1 ? 'pixel' : 'instance';

                const label = e.label || e.file.split('.')[0];
                const sizeFactor = e.sizeFactor || 1;

                build = build
                    .toRoot()
                    .apply(ReadFile, { file, label })
                    .apply(TrajectoryFromGRO)
                    .apply(ModelFromTrajectory, { modelIndex: 0 })
                    .apply(StructureFromGeneric, { transforms, label })
                    .apply(StructureRepresentation3D, getSpacefillParams(color, sizeFactor, customState.lodLevels, clipVariant), { tags: e.groups });
            }
            await build.commit();

            const rootId = `${manifest.roots[0].id}:`;
            const values = { type: 'group-generate', value: ColorNames.white, lightness: 0, alpha: 1 };
            await updateColors(ctx, values, rootId, '');
        } catch (e) {
            console.error(e);
            ctx.log.error(e);
        } finally {
            ctx.animationLoop.start();
        }
    }).run();
}

//

type GenericRoot = {
    id: string
    label?: string
    description?: string
    children: GenericGroup[]
}

type GenericGroup = {
    id: string
    root: string
    label?: string
    description?: string
    children?: GenericGroup[]
}

type GenericEntity = {
    file: string
    label?: string
    description?: string
    /** references to `${GenericGroup.id}:${GenericGroup.root}` */
    groups: string[]
    instances: GenericInstances
    /** defaults to 1 */
    sizeFactor?: number
}

type GenericInstances = {
    /** translation vectors [x0, y0, z0, ..., xn, yn, zn] with n = count - 1 */
    positions: number[]
    /** euler angles [x0, y0, z0, ..., xn, yn, zn] with n = count - 1 */
    rotations: number[]
}

type GenericManifest = {
    label?: string
    description?: string
    roots: GenericRoot[]
    entities: GenericEntity[]
}

//

type MartiniManifest = {
    model: string,
    positions: Vec3[],
    rotations: Vec3[],
    function: string
}[]

function martiniToGeneric(martini: MartiniManifest): GenericManifest {
    const functionRoot: GenericRoot = {
        id: 'function',
        label: 'Function',
        description: 'Functional classification',
        children: [],
    };
    const entities: GenericEntity[] = [];

    const seenGroups = new Set<string>();

    const membraneGroup = {
        id: 'membane',
        root: 'function',
        label: 'Membrane',
        children: [] as GenericGroup[],
    };
    functionRoot.children!.push(membraneGroup);
    seenGroups.add(membraneGroup.id);

    const lipidsGroup = {
        id: 'lipid',
        root: 'function',
        label: 'Lipid',
        children: [] as GenericGroup[],
    };
    membraneGroup.children!.push(lipidsGroup);
    seenGroups.add(lipidsGroup.id);

    const upperGroup = {
        id: 'upper',
        root: 'function',
        label: 'Upper Leaflet',
    };
    lipidsGroup.children!.push(upperGroup);
    seenGroups.add(upperGroup.id);

    const lowerGroup = {
        id: 'lower',
        root: 'function',
        label: 'Lower Leaflet',
    };
    lipidsGroup.children!.push(lowerGroup);
    seenGroups.add(lowerGroup.id);

    const memprotGroup = {
        id: 'memprot',
        root: 'function',
        label: 'Transmembrane Protein',
    };
    membraneGroup.children!.push(memprotGroup);
    seenGroups.add(memprotGroup.id);

    for (const e of martini) {
        const label = e.model.split('.')[0];
        const group = e.function || 'Metabolite';

        const positions = e.positions.flat().map(x => Math.round((x * 10) * 100) / 100);
        const rotations = e.rotations.flat().map(x => Math.round(x * 100) / 100);

        if (group.includes('lower leaflet')) {
            entities.push({
                file: e.model,
                label: label.substring(15),
                groups: [`function:lower`],
                instances: { positions, rotations },
                sizeFactor: 1.5,
            });
        } else if (group.includes('upper leaflet')) {
            entities.push({
                file: e.model,
                label: label.substring(15),
                groups: [`function:upper`],
                instances: { positions, rotations },
                sizeFactor: 1.5,
            });
        } else if (group.length === 4) {
            entities.push({
                file: e.model,
                label: label.substring(17),
                groups: [`function:memprot`],
                instances: { positions, rotations },
                sizeFactor: 1.5,
            });
        } else {
            if (!seenGroups.has(group)) {
                functionRoot.children!.push({
                    id: group,
                    root: 'function',
                    label: group,
                });
                seenGroups.add(group);
            }
            entities.push({
                file: e.model,
                label,
                groups: [`function:${group}`],
                instances: { positions, rotations },
                sizeFactor: 1.5,
            });
        }
    }

    return {
        label: 'Martini',
        description: 'Martini coarse-grained model',
        roots: [functionRoot],
        entities,
    };
}
