/**
 * Copyright (c) 2023-2024 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4 } from '../../../../mol-math/linear-algebra/3d/mat4';
import { StateBuilder, StateObjectSelector } from '../../../../mol-state';
import { PluginContext } from '../../../../mol-plugin/context';
import { SpacefillRepresentationProvider } from '../../../../mol-repr/structure/representation/spacefill';
import { Color } from '../../../../mol-util/color';
import { utf8Read } from '../../../../mol-io/common/utf8';
import { Mat3, Quat, Vec3 } from '../../../../mol-math/linear-algebra';
import { GraphicsMode, MesoscaleGroup, MesoscaleState, getGraphicsModeProps, getMesoscaleGroupParams } from '../state';
import { ColorNames } from '../../../../mol-util/color/names';
import { ShapeRepresentation3D, StructureRepresentation3D } from '../../../../mol-plugin-state/transforms/representation';
import { ParseCif, ParsePly, ReadFile } from '../../../../mol-plugin-state/transforms/data';
import { ModelFromTrajectory, ShapeFromPly, TrajectoryFromGRO, TrajectoryFromMOL, TrajectoryFromMOL2, TrajectoryFromMmCif, TrajectoryFromPDB, TrajectoryFromSDF, TrajectoryFromXYZ } from '../../../../mol-plugin-state/transforms/model';
import { Euler } from '../../../../mol-math/linear-algebra/3d/euler';
import { Asset } from '../../../../mol-util/assets';
import { Clip } from '../../../../mol-util/clip';
import { StructureFromGeneric } from './model';
import { getFileNameInfo } from '../../../../mol-util/file-info';
import { NumberArray } from '../../../../mol-util/type-helpers';
import { BaseGeometry } from '../../../../mol-geo/geometry/base';
import { ParamDefinition as PD } from '../../../../mol-util/param-definition';

function getSpacefillParams(color: Color, sizeFactor: number, graphics: GraphicsMode, clipVariant: Clip.Variant) {
    const gmp = getGraphicsModeProps(graphics === 'custom' ? 'quality' : graphics);
    return {
        type: {
            name: 'spacefill',
            params: {
                ...SpacefillRepresentationProvider.defaultValues,
                ignoreHydrogens: true,
                instanceGranularity: true,
                ignoreLight: true,
                lodLevels: gmp.lodLevels.map(l => {
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

function getPlyShapeParams(color: Color, clipVariant: Clip.Variant) {
    return {
        ...PD.getDefaultValues(BaseGeometry.Params),
        instanceGranularity: true,
        ignoreLight: true,
        clip: {
            variant: clipVariant,
            objects: [],
        },
        quality: 'custom',
        doubleSided: true,
        coloring: {
            name: 'uniform',
            params: { color }
        },
        grouping: {
            name: 'none',
            params: {}
        },
        material: {
            metalness: 0.0,
            roughness: 1.0,
            bumpiness: 1.0,
        },
        bumpAmplitude: 0.1,
        bumpFrequency: 0.1 / 10,
    };
}

export async function createGenericHierarchy(plugin: PluginContext, file: Asset.File) {
    const asset = await plugin.runTask(plugin.managers.asset.resolve(file, 'zip'));
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

    const state = plugin.state.data;
    const graphicsMode = MesoscaleState.get(plugin).graphics;
    const groupParams = getMesoscaleGroupParams(graphicsMode);

    async function addGroup(g: GenericGroup, cell: StateObjectSelector, parent: string) {
        const group = await state.build()
            .to(cell)
            .apply(MesoscaleGroup, { ...groupParams, index: undefined, tag: `${g.root}:${g.id}`, label: g.label || g.id, description: g.description, color: { type: 'custom', illustrative: false, value: ColorNames.white, variability: 20, shift: 0, lightness: 0, alpha: 1, emissive: 0 } }, { tags: [`group:${g.root}:${g.id}`, g.root === parent ? `${g.root}:` : `${g.root}:${parent}`], state: { isCollapsed: true, isHidden: groupParams.hidden } })
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
            .apply(MesoscaleGroup, { ...groupParams, root: true, index: -1, tag: `${r.id}:`, label: r.label || r.id, description: r.description, color: { type: 'custom', illustrative: false, value: ColorNames.white, variability: 20, shift: 0, lightness: 0, alpha: 1, emissive: 0 } }, { tags: `group:${r.id}:`, state: { isCollapsed: false, isHidden: groupParams.hidden } })
            .commit();

        if (r.children) {
            for (const c of r.children!) {
                await addGroup(c, root, r.id);
            }
        }
    }

    const transformAssets = new Map<string, Asset>();
    const getTransformAsset = (file: string) => {
        if (!transformAssets.has(file)) {
            const d = asset.data[file];
            transformAssets.set(file, Asset.File(new File([d], file)));
        }
        return transformAssets.get(file)!;
    };

    const getAssetInstances = (instances: GenericInstances<string>): GenericInstances<Asset> => {
        return {
            positions: {
                data: Array.isArray(instances.positions.data)
                    ? instances.positions.data
                    : {
                        file: getTransformAsset(instances.positions.data.file),
                        view: instances.positions.data.view,
                    },
                type: instances.positions.type,
            },
            rotations: {
                data: Array.isArray(instances.rotations.data)
                    ? instances.rotations.data
                    : {
                        file: getTransformAsset(instances.rotations.data.file),
                        view: instances.rotations.data.view,
                    },
                variant: instances.rotations.variant,
                type: instances.rotations.type,
            }
        };
    };

    await state.transaction(async () => {
        try {
            plugin.animationLoop.stop({ noDraw: true });
            let build: StateBuilder.Root | StateBuilder.To<any> = state.build();
            for (const ent of manifest.entities) {
                const d = asset.data[ent.file];
                const info = getFileNameInfo(ent.file);
                const isBinary = ['bcif'].includes(info.ext);

                const t = isBinary ? d : utf8Read(d, 0, d.length);
                const file = Asset.File(new File([t], ent.file));

                const color = ColorNames.skyblue;

                const sizeFactor = ent.sizeFactor || 1;
                const tags = ent.groups.map(({ id, root }) => `${root}:${id}`);
                const instances = ent.instances && getAssetInstances(ent.instances);
                const description = ent.description;
                const label = ent.label || ent.file.split('.')[0];
                build = build
                    .toRoot()
                    .apply(ReadFile, { file, label, isBinary });

                if (['gro', 'cif', 'mmcif', 'mcif', 'bcif', 'pdb', 'ent', 'xyz', 'mol', 'sdf', 'sd', 'mol2'].includes(info.ext)) {
                    if (['gro'].includes(info.ext)) {
                        build = build.apply(TrajectoryFromGRO);
                    } else if (['cif', 'mmcif', 'mcif', 'bcif'].includes(info.ext)) {
                        build = build.apply(ParseCif).apply(TrajectoryFromMmCif);
                    } else if (['pdb', 'ent'].includes(info.ext)) {
                        build = build.apply(TrajectoryFromPDB);
                    } else if (['xyz'].includes(info.ext)) {
                        build = build.apply(TrajectoryFromXYZ);
                    } else if (['mol'].includes(info.ext)) {
                        build = build.apply(TrajectoryFromMOL);
                    } else if (['sdf', 'sd'].includes(info.ext)) {
                        build = build.apply(TrajectoryFromSDF);
                    } else if (['mol2'].includes(info.ext)) {
                        build = build.apply(TrajectoryFromMOL2);
                    }

                    let clipVariant: Clip.Variant = 'pixel';
                    if (ent.instances) {
                        if (Array.isArray(ent.instances.positions.data)) {
                            clipVariant = ent.instances.positions.data.length <= 3 ? 'pixel' : 'instance';
                        } else {
                            const byteLength = ent.instances.positions.data.view
                                ? ent.instances.positions.data.view.byteLength
                                : asset.data[ent.instances.positions.data.file].length;
                            clipVariant = byteLength <= 3 * 4 ? 'pixel' : 'instance';
                        }
                    }

                    build = build
                        .apply(ModelFromTrajectory, { modelIndex: 0 })
                        .apply(StructureFromGeneric, { instances, label, description })
                        .apply(StructureRepresentation3D, getSpacefillParams(color, sizeFactor, graphicsMode, clipVariant), { tags });
                } else if (['ply'].includes(info.ext)) {
                    if (['ply'].includes(info.ext)) {
                        const transforms = await getTransforms(plugin, instances);
                        const clipVariant = transforms.length === 1 ? 'pixel' : 'instance';
                        build = build
                            .apply(ParsePly)
                            .apply(ShapeFromPly, { label, transforms })
                            .apply(ShapeRepresentation3D, getPlyShapeParams(color, clipVariant), { tags });
                    }
                } else {
                    console.warn(`unknown file format '${info.ext}'`);
                }
            }
            await build.commit();
        } catch (e) {
            console.error(e);
            plugin.log.error(e);
        } finally {
            plugin.animationLoop.start();
        }
    }).run();

    asset.dispose();
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
    /** reference to `${GenericRoot.id}` */
    root: string
    label?: string
    description?: string
    children?: GenericGroup[]
}

type GenericEntity = {
    /**
     * the entity file name
     *
     * the following extensions/formats are supported
     *
     * structures
     * - gro
     * - cif, mmcif, mcif, bcif
     * - pdb, ent
     * - xyz
     * - mol
     * - sdf, sd
     * - mol2
     *
     * meshes
     * - ply
     */
    file: string
    label?: string
    description?: string
    groups: {
        /** reference to `${GenericGroup.id}` */
        id: string,
        /** reference to `${GenericGroup.root}` */
        root: string
    }[]
    /**
     * defaults to a single, untransformed instance
     */
    instances?: GenericInstances<string>
    /**
     * defaults to 1 (assuming fully atomic structures)
     * for C-alpha only structures set to 2
     * for Martini coarse-grained set to 1.5
     */
    sizeFactor?: number
}

type BinaryData<T extends string | Asset> = {
    file: T,
    view?: {
        byteOffset: number,
        byteLength: number
    }
}

export type GenericInstances<T extends string | Asset> = {
    /**
     * translation vectors in Angstrom
     * [x0, y0, z0, ..., xn, yn, zn] with n = count - 1
     */
    positions: {
        /**
         * either the data itself or a pointer to binary data
         */
        data: number[] | BinaryData<T>
        /**
         * how to interpret the data
         * defaults to `{ kind: 'Array', type: 'Float32' }`
         */
        type?: { kind: 'Array', type: 'Float32' }
        // TODO: maybe worthwhile in the future, mirroring encoders from BinaryCIF
        // | { kind: 'IntegerPackedFixedPoint', byteCount: number, srcSize: number, factor: number, srcType: 'Float32' }
    }
    /**
     * euler angles in XYZ order
     * [x0, y0, z0, ..., xn, yn, zn] with n = count - 1
     *
     * quaternion rotations in XYZW order
     * [x0, y0, z0, w0, ..., xn, yn, zn, wn] with n = count - 1
     *
     * rotation matrices in row-major order
     * [m00_0, m01_0, m02_0, ..., m20_n, m21_n, m22_n] with n = count - 1
     */
    rotations: {
        variant: 'euler' | 'quaternion' | 'matrix',
        /**
         * either the data itself or a pointer to binary data
         */
        data: number[] | BinaryData<T>
        /**
         * how to interpret the data
         * defaults to `{ kind: 'Array', type: 'Float32' }`
         */
        type?: { kind: 'Array', type: 'Float32' }
    }
}

type GenericFrame = {
    time: number
    entities: {
        file: string
        instances: GenericInstances<string>
    }[]
}

type GenericTrajectory = {
    label?: string
    description?: string
    frames: GenericFrame[]
}

type GenericManifest = {
    label?: string
    description?: string
    roots: GenericRoot[]
    entities: GenericEntity[]
    trajectories?: GenericTrajectory[]
}

//

const p = Vec3();
const q = Quat();
const m = Mat3();
const e = Euler();

async function getPositions(plugin: PluginContext, p: GenericInstances<Asset>['positions']): Promise<NumberArray> {
    if (Array.isArray(p.data)) {
        return p.data;
    } else {
        const a = await plugin.runTask(plugin.managers.asset.resolve(p.data.file, 'binary'));
        const o = p.data.view?.byteOffset || 0;
        const l = p.data.view?.byteLength || a.data.byteLength;
        return new Float32Array(a.data.buffer, o + a.data.byteOffset, l / 4);
    }
};

async function getRotations(plugin: PluginContext, r: GenericInstances<Asset>['rotations']): Promise<NumberArray> {
    if (Array.isArray(r.data)) {
        return r.data;
    } else {
        const a = await plugin.runTask(plugin.managers.asset.resolve(r.data.file, 'binary'));
        const o = r.data.view?.byteOffset || 0;
        const l = r.data.view?.byteLength || a.data.byteLength;
        return new Float32Array(a.data.buffer, o + a.data.byteOffset, l / 4);
    }
};

export async function getTransforms(plugin: PluginContext, instances?: GenericInstances<Asset>) {
    const transforms: Mat4[] = [];
    if (instances) {
        const positions = await getPositions(plugin, instances.positions);
        const rotations = await getRotations(plugin, instances.rotations);

        for (let i = 0, il = positions.length / 3; i < il; ++i) {
            Vec3.fromArray(p, positions, i * 3);
            if (instances.rotations.variant === 'matrix') {
                Mat3.fromArray(m, rotations, i * 9);
                const t = Mat4.fromMat3(Mat4(), m);
                Mat4.setTranslation(t, p);
                transforms.push(t);
            } else if (instances.rotations.variant === 'quaternion') {
                Quat.fromArray(q, rotations, i * 4);
                const t = Mat4.fromQuat(Mat4(), q);
                Mat4.setTranslation(t, p);
                transforms.push(t);
            } else if (instances.rotations.variant === 'euler') {
                Euler.fromArray(e, rotations, i * 3);
                Quat.fromEuler(q, e, 'XYZ');
                const t = Mat4.fromQuat(Mat4(), q);
                Mat4.setTranslation(t, p);
                transforms.push(t);
            }
        }
    } else {
        transforms.push(Mat4.identity());
    }

    return transforms;
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

        const positions = {
            data: e.positions.flat().map(x => Math.round((x * 10) * 100) / 100)
        };
        const rotations = {
            data: e.rotations.flat().map(x => Math.round(x * 100) / 100),
            variant: 'euler' as const,
        };

        if (group.includes('lower leaflet')) {
            entities.push({
                file: e.model,
                label: label.substring(15),
                groups: [{ root: 'function', id: 'lower' }],
                instances: { positions, rotations },
                sizeFactor: 1.5,
            });
        } else if (group.includes('upper leaflet')) {
            entities.push({
                file: e.model,
                label: label.substring(15),
                groups: [{ root: 'function', id: 'upper' }],
                instances: { positions, rotations },
                sizeFactor: 1.5,
            });
        } else if (group.length === 4) {
            entities.push({
                file: e.model,
                label: label.substring(17),
                groups: [{ root: 'function', id: 'memprot' }],
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
                groups: [{ root: 'function', id: group }],
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
