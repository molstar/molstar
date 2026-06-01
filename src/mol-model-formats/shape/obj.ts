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
import { deepClone } from '../../mol-util/object';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { distinctColors } from '../../mol-util/color/distinct';
import { ValueCell } from '../../mol-util/value-cell';

export type ObjData = {
    source: ObjFile,
    transforms?: Mat4[],
}

export const ObjShapeParams = {
    ...Mesh.Params,
    coloring: PD.MappedStatic('group', {
        group: PD.Group({}),
        uniform: PD.Group({
            color: PD.Color(ColorNames.grey),
        }, { isFlat: true })
    }),
};
export type ObjShapeParams = typeof ObjShapeParams

/**
 * Build an expanded mesh from indexed OBJ data.
 *
 * OBJ stores positions and normals as indexed arrays (each face-vertex
 * references a position index and an optional normal index). Mesh.create
 * requires a flat vertex array, so we expand the indexed data.
 *
 * A vertex key encodes (posIdx, normIdx, groupIdx) so that a position
 * shared between faces of different material groups or with different normals
 * gets a distinct mesh vertex — ensuring each vertex has one unambiguous
 * group ID and normal.
 */
async function getMesh(ctx: RuntimeContext, obj: ObjFile, mesh?: Mesh): Promise<Mesh> {
    const { positions, normals, positionIndices, normalIndices, groupIndices, triangleCount } = obj;
    const hasNormals = obj.normalCount > 0;

    const builderState = MeshBuilder.createState(triangleCount * 3, triangleCount, mesh);
    const { vertices, normals: normBuf, indices, groups } = builderState;

    // Group count for the composite dedup key; groups always contains at least 'default'.
    const groupCount = Math.max(1, obj.groups.length);

    // Map from position index to a map of (normal index, group) composite keys to expanded vertex index.
    // Avoids allocating a string key per face-vertex on the hot path.
    const vertexMap = new Map<number, Map<number, number>>();

    const getVertex = (pi: number, ni: number, groupId: number): number => {
        let inner = vertexMap.get(pi);
        if (inner === undefined) {
            inner = new Map<number, number>();
            vertexMap.set(pi, inner);
        }
        const subKey = (ni + 1) * groupCount + groupId;
        let idx = inner.get(subKey);
        if (idx === undefined) {
            idx = vertices.elementCount;
            inner.set(subKey, idx);

            const po = pi * 3;
            ChunkedArray.add3(vertices, positions[po], positions[po + 1], positions[po + 2]);

            if (hasNormals && ni >= 0) {
                const no = ni * 3;
                ChunkedArray.add3(normBuf, normals[no], normals[no + 1], normals[no + 2]);
            } else {
                ChunkedArray.add3(normBuf, 0, 0, 0);
            }

            ChunkedArray.add(groups, groupId);
        }
        return idx;
    };

    const updateChunk = 50000;

    for (let t = 0; t < triangleCount; ++t) {
        const triOffset = t * 3;
        const groupId = groupIndices[t];

        const i0 = getVertex(positionIndices[triOffset], hasNormals ? normalIndices[triOffset] : -1, groupId);
        const i1 = getVertex(positionIndices[triOffset + 1], hasNormals ? normalIndices[triOffset + 1] : -1, groupId);
        const i2 = getVertex(positionIndices[triOffset + 2], hasNormals ? normalIndices[triOffset + 2] : -1, groupId);

        ChunkedArray.add3(indices, i0, i1, i2);

        if (t % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Building OBJ mesh', current: t, max: triangleCount });
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    if (!hasNormals) Mesh.computeNormals(m);

    ValueCell.updateIfChanged(m.varyingGroup, true);

    return m;
}

type GroupColors = { kind: 'group', colors: Color[] } | { kind: 'uniform', color: Color }

function getColoring(obj: ObjFile, props: PD.Values<ObjShapeParams>): GroupColors {
    const { coloring } = props;
    if (coloring.name === 'uniform') {
        return { kind: 'uniform', color: coloring.params.color };
    }
    // Generate distinct colors for each group
    const n = Math.max(1, obj.groups.length);
    const colors = distinctColors(n);
    return { kind: 'group', colors };
}

function createShape(objData: ObjData, mesh: Mesh, groupColors: GroupColors) {
    const { source, transforms } = objData;
    return Shape.create(
        'obj-mesh', source, mesh,
        (groupId: number) => {
            if (groupColors.kind === 'uniform') return groupColors.color;
            const idx = Math.min(groupId, groupColors.colors.length - 1);
            return groupColors.colors[idx];
        },
        () => 1,
        (groupId: number) => {
            const name = source.groups[groupId] ?? `Group ${groupId}`;
            return name;
        },
        transforms
    );
}

function makeShapeGetter() {
    let _objData: ObjData | undefined;
    let _props: PD.Values<ObjShapeParams> | undefined;

    let _shape: Shape<Mesh>;
    let _mesh: Mesh;
    let _groupColors: GroupColors;

    const getShape = async (ctx: RuntimeContext, objData: ObjData, props: PD.Values<ObjShapeParams>, shape?: Shape<Mesh>) => {
        let newMesh = false;
        let newColor = false;

        if (!_objData || _objData !== objData) {
            newMesh = true;
        }

        if (!_props || !PD.isParamEqual(ObjShapeParams.coloring, _props.coloring, props.coloring)) {
            newColor = true;
        }

        if (newMesh) {
            _mesh = await getMesh(ctx, objData.source, shape && shape.geometry);
            _groupColors = getColoring(objData.source, props);
            _shape = createShape(objData, _mesh, _groupColors);
        } else if (newColor) {
            _groupColors = getColoring(objData.source, props);
            _shape = createShape(objData, _mesh, _groupColors);
        }

        _objData = objData;
        _props = deepClone(props);

        return _shape;
    };
    return getShape;
}

export function shapeFromObj(source: ObjFile, params?: { transforms?: Mat4[] }) {
    return Task.create<ShapeProvider<ObjData, Mesh, ObjShapeParams>>('Shape Provider', async _ctx => {
        return {
            label: 'Mesh',
            data: { source, transforms: params?.transforms },
            params: ObjShapeParams,
            getShape: makeShapeGetter(),
            geometryUtils: Mesh.Utils
        };
    });
}
