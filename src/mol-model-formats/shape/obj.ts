/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { ObjFile } from '../../mol-io/reader/obj/schema';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Shape } from '../../mol-model/shape';
import { ChunkedArray } from '../../mol-data/util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { distinctColors } from '../../mol-util/color/distinct';

export type ObjData = {
    source: ObjFile,
    transforms?: Mat4[],
}

function createObjShapeParams(objFile?: ObjFile) {
    const materialNames = objFile?.materialNames ?? [];
    const hasMaterials = materialNames.length > 0;

    const defaultColors = materialNames.length > 1
        ? distinctColors(materialNames.length)
        : materialNames.length === 1 ? [ColorNames.grey] : [];

    const materialColorParams: Record<string, PD.Color> = {};
    for (let i = 0; i < materialNames.length; ++i) {
        materialColorParams[materialNames[i]] = PD.Color(defaultColors[i]);
    }

    return {
        ...Mesh.Params,
        coloring: PD.MappedStatic(hasMaterials ? 'material' : 'uniform', {
            uniform: PD.Group({
                color: PD.Color(ColorNames.grey),
            }, { isFlat: true }),
            material: PD.Group(materialColorParams, { isFlat: false }),
        }),
    };
}

export const ObjShapeParams = createObjShapeParams();
export type ObjShapeParams = typeof ObjShapeParams

/**
 * Resolve one Color per material group from the current params.
 * The returned array has length = max(1, materialNames.length).
 */
function getMaterialColors(materialNames: readonly string[], props: PD.Values<ObjShapeParams>): Color[] {
    const count = Math.max(1, materialNames.length);
    const { coloring } = props;
    if (coloring.name === 'uniform') {
        return Array<Color>(count).fill(coloring.params.color);
    } else {
        // material: read one color per material name from dynamic params
        if (materialNames.length === 0) return [ColorNames.grey];
        const params = coloring.params as Record<string, Color>;
        return materialNames.map(name => params[name] ?? ColorNames.grey);
    }
}

/**
 * Build a mesh from indexed OBJ data.
 *
 * OBJ stores positions and normals as indexed arrays (each face-vertex
 * references a position index and an optional normal index). Mesh.create
 * requires a flat vertex array, so we expand the indexed data.
 *
 * A vertex key encodes (posIdx, normIdx) so that a shared position with
 * different normals gets a distinct mesh vertex.
 */
async function getMesh(ctx: RuntimeContext, obj: ObjFile, mesh?: Mesh): Promise<Mesh> {
    const { positions, normals, positionIndices, normalIndices, triangleCount, faceGroups } = obj;
    const hasNormals = obj.normalCount > 0;

    const builderState = MeshBuilder.createState(triangleCount * 3, triangleCount, mesh);
    const { vertices, normals: normBuf, indices, groups } = builderState;

    const updateChunk = 50000;

    for (let t = 0; t < triangleCount; ++t) {
        const triOffset = t * 3;
        const base = t * 3;

        for (let v = 0; v < 3; ++v) {
            const pi = positionIndices[triOffset + v];
            const po = pi * 3;
            ChunkedArray.add3(vertices, positions[po], positions[po + 1], positions[po + 2]);

            const ni = hasNormals ? normalIndices[triOffset + v] : -1;
            if (hasNormals && ni >= 0) {
                const no = ni * 3;
                ChunkedArray.add3(normBuf, normals[no], normals[no + 1], normals[no + 2]);
            } else {
                ChunkedArray.add3(normBuf, 0, 0, 0);
            }

            ChunkedArray.add(groups, faceGroups[t]);
        }

        ChunkedArray.add3(indices, base, base + 1, base + 2);

        if (t % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Building OBJ mesh', current: t, max: triangleCount });
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    if (!hasNormals) Mesh.computeNormals(m);

    return m;
}

function createShape(objData: ObjData, mesh: Mesh, colors: Color[]) {
    const { source, transforms } = objData;
    const { materialNames } = source;
    return Shape.create(
        'obj-mesh', source, mesh,
        (groupId: number) => colors[Math.min(groupId, colors.length - 1)],
        () => 1,
        (groupId: number) => materialNames.length > 0 ? (materialNames[groupId] ?? 'OBJ Mesh') : 'OBJ Mesh',
        transforms,
        colors.length
    );
}

function makeShapeGetter() {
    let _objData: ObjData | undefined;
    let _colors: Color[] | undefined;

    let _shape: Shape<Mesh>;
    let _mesh: Mesh;

    const getShape = async (ctx: RuntimeContext, objData: ObjData, props: PD.Values<ObjShapeParams>, shape?: Shape<Mesh>) => {
        const newMesh = !_objData || _objData !== objData;

        const nextColors = getMaterialColors(objData.source.materialNames, props);
        const newColor = !_colors
            || nextColors.length !== _colors.length
            || nextColors.some((c, i) => c !== _colors![i]);

        if (newMesh) {
            _colors = nextColors;
            _mesh = await getMesh(ctx, objData.source, shape && shape.geometry);
            _shape = createShape(objData, _mesh, _colors);
        } else if (newColor) {
            _colors = nextColors;
            _shape = createShape(objData, _mesh, _colors);
        }

        _objData = objData;

        return _shape;
    };
    return getShape;
}

export function shapeFromObj(source: ObjFile, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<ObjData, Mesh, ObjShapeParams>>('Shape Provider', async _ctx => {
        return {
            label: 'Mesh',
            data: { source, transforms: params?.transforms },
            params: createObjShapeParams(source),
            getShape: makeShapeGetter(),
            geometryUtils: Mesh.Utils,
        };
    });
}

