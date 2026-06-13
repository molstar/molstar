/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext, Task } from '../../mol-task';
import { ShapeProvider } from '../../mol-model/shape/provider';
import { Color } from '../../mol-util/color';
import { MtlFile, ObjFile } from '../../mol-io/reader/obj/schema';
import { MeshBuilder } from '../../mol-geo/geometry/mesh/mesh-builder';
import { Mesh } from '../../mol-geo/geometry/mesh/mesh';
import { Shape } from '../../mol-model/shape';
import { ChunkedArray } from '../../mol-data/util';
import { ParamDefinition as PD } from '../../mol-util/param-definition';
import { ColorNames } from '../../mol-util/color/names';
import { Mat4 } from '../../mol-math/linear-algebra/3d/mat4';
import { distinctColors } from '../../mol-util/color/distinct';
import { ValueCell } from '../../mol-util';
import { isDebugMode } from '../../mol-util/debug';

export type ObjData = {
    source: ObjFile,
    transforms?: Mat4[],
    mtl?: MtlFile,
}

type ColoringName = 'uniform' | 'given' | 'custom' | 'vertex';
type ColoringParams = {
    uniform: PD.Group<{ color: Color }>,
    given?: PD.Group<{}>,
    custom?: PD.Group<Record<string, Color>>,
    vertex?: PD.Group<{}>,
}

export function createObjShapeParams(objFile?: ObjFile, mtl?: MtlFile) {
    const hasVertexColors = (objFile?.vertexColors.length ?? 0) > 0;
    const materialNames = objFile?.materialNames ?? [];
    const hasMaterials = materialNames.length > 0;
    const hasMtl = mtl !== undefined && mtl.size > 0 && hasMaterials;

    const defaultColors = materialNames.length > 1
        ? distinctColors(materialNames.length)
        : materialNames.length === 1 ? [ColorNames.grey] : [];

    const materialColorParams: Record<string, PD.Color> = {};
    for (let i = 0; i < materialNames.length; ++i) {
        materialColorParams[materialNames[i]] = PD.Color(defaultColors[i]);
    }

    const coloringOptions: ColoringParams = {
        uniform: PD.Group({
            color: PD.Color(ColorNames.grey),
        }, { isFlat: true }),
    };
    if (hasMtl) {
        coloringOptions.given = PD.Group({}, { isFlat: true });
    }
    coloringOptions.custom = PD.Group(materialColorParams, { isFlat: false });
    if (hasVertexColors) {
        coloringOptions.vertex = PD.Group({}, { isFlat: true });
    }
    const defaultColoring: ColoringName = hasVertexColors ? 'vertex' : hasMtl ? 'given' : hasMaterials ? 'custom' : 'uniform';

    return {
        ...Mesh.Params,
        coloring: PD.MappedStatic(defaultColoring, coloringOptions),
    };
}

export const ObjShapeParams = createObjShapeParams();
export type ObjShapeParams = typeof ObjShapeParams

/**
 * Resolve one Color per material group from the current params.
 * Returns an empty array when vertex-color mode is active (colors come from geometry).
 */
function getMaterialColors(materialNames: readonly string[], props: PD.Values<ObjShapeParams>, mtl?: MtlFile): Color[] {
    const { coloring } = props;
    if (coloring.name === 'vertex') return [];
    const count = Math.max(1, materialNames.length);
    if (coloring.name === 'uniform') {
        return Array<Color>(count).fill(coloring.params.color);
    } else if (coloring.name === 'given') {
        // given: read Kd directly from the MTL file
        if (materialNames.length === 0) return [ColorNames.grey];
        return materialNames.map(name => mtl?.get(name)?.Kd ?? ColorNames.grey);
    } else {
        // custom: read one color per material name from user-editable params
        if (materialNames.length === 0) return [ColorNames.grey];
        const params = coloring.params as Record<string, Color>;
        return materialNames.map(name => params[name] ?? ColorNames.grey);
    }
}

/**
 * Build a per-position → material-index lookup for use when the mesh is built
 * with position-index groups (i.e. when vertex colors are present).
 */
function buildPositionToMaterial(obj: ObjFile): Int32Array {
    const map = new Int32Array(obj.positionCount);
    for (let t = 0; t < obj.triangleCount; ++t) {
        const triOffset = t * 3;
        const g = obj.faceGroups[t];
        for (let v = 0; v < 3; ++v) {
            map[obj.positionIndices[triOffset + v]] = g;
        }
    }
    return map;
}

function getCanUseIndexedPath(obj: ObjFile, usePositionGroups: boolean): { canUseIndexedPath: boolean, posGroups: Int32Array | null, posNormals: Int32Array | null } {
    const { positionIndices, normalIndices, faceGroups, positionCount, triangleCount } = obj;
    const hasNormals = obj.normalCount > 0;
    let canUseIndexedPath = usePositionGroups;
    let posGroups: Int32Array | null = null;
    let posNormals: Int32Array | null = null;

    if (!usePositionGroups) {
        posGroups = new Int32Array(positionCount).fill(-1);
        let hasGroupConflict = false;
        outer: for (let t = 0; t < triangleCount; ++t) {
            const triOffset = t * 3;
            const g = faceGroups[t];
            for (let v = 0; v < 3; ++v) {
                const pi = positionIndices[triOffset + v];
                if (posGroups[pi] === -1) {
                    posGroups[pi] = g;
                } else if (posGroups[pi] !== g) {
                    hasGroupConflict = true;
                    break outer;
                }
            }
        }
        if (hasGroupConflict) {
            if (isDebugMode) console.warn('OBJ: position shared across face groups; falling back to vertex duplication');
        } else {
            canUseIndexedPath = true;
        }
    }

    if (canUseIndexedPath && hasNormals) {
        posNormals = new Int32Array(positionCount).fill(-1);
        let hasNormalConflict = false;
        outerNormals: for (let t = 0; t < triangleCount; ++t) {
            const triOffset = t * 3;
            for (let v = 0; v < 3; ++v) {
                const pi = positionIndices[triOffset + v];
                const ni = normalIndices[triOffset + v];
                if (posNormals[pi] === -1) {
                    posNormals[pi] = ni;
                } else if (posNormals[pi] !== ni) {
                    hasNormalConflict = true;
                    break outerNormals;
                }
            }
        }
        if (hasNormalConflict) {
            canUseIndexedPath = false;
        }
    }
    return { canUseIndexedPath, posGroups, posNormals };
}

type IndexedInput = {
    posGroups: Int32Array | null
    posNormals: Int32Array | null
    usePositionGroups: boolean
}

async function getIndexedMesh(ctx: RuntimeContext, obj: ObjFile, input: IndexedInput, mesh?: Mesh): Promise<Mesh> {
    const { positions, normals, positionIndices, triangleCount, positionCount } = obj;
    const hasNormals = obj.normalCount > 0;

    const { posGroups, posNormals, usePositionGroups } = input;

    const updateChunk = 50000;

    const builderState = MeshBuilder.createState(positionCount, triangleCount, mesh);
    const { vertices, normals: normBuf, indices, groups } = builderState;

    for (let p = 0; p < positionCount; ++p) {
        const po = p * 3;
        ChunkedArray.add3(vertices, positions[po], positions[po + 1], positions[po + 2]);

        if (posNormals !== null && posNormals[p] >= 0) {
            const no = posNormals[p] * 3;
            ChunkedArray.add3(normBuf, normals[no], normals[no + 1], normals[no + 2]);
        } else {
            ChunkedArray.add3(normBuf, 0, 0, 0);
        }

        // position-index groups: group = position index; material groups: group = material idx
        ChunkedArray.add(groups, usePositionGroups ? p : Math.max(0, posGroups![p]));

        if (p % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Building OBJ mesh', current: p, max: positionCount });
        }
    }

    for (let t = 0; t < triangleCount; ++t) {
        const triOffset = t * 3;
        ChunkedArray.add3(indices,
            positionIndices[triOffset],
            positionIndices[triOffset + 1],
            positionIndices[triOffset + 2]
        );
    }

    const m = MeshBuilder.getMesh(builderState);
    if (!hasNormals) Mesh.computeNormals(m);

    // TODO: check if needed
    ValueCell.updateIfChanged(m.varyingGroup, true);

    return m;
}

async function getFaceMesh(ctx: RuntimeContext, obj: ObjFile, usePositionGroups: boolean, mesh?: Mesh): Promise<Mesh> {
    const { positions, normals, positionIndices, normalIndices, triangleCount, faceGroups } = obj;
    const hasNormals = obj.normalCount > 0;

    const updateChunk = 50000;

    const builderState = MeshBuilder.createState(triangleCount * 3, triangleCount, mesh);
    const { vertices, normals: normBuf, indices, groups } = builderState;

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

            // position-index groups: group = position index; material groups: group = material idx
            ChunkedArray.add(groups, usePositionGroups ? pi : faceGroups[t]);
        }

        ChunkedArray.add3(indices, base, base + 1, base + 2);

        if (t % updateChunk === 0 && ctx.shouldUpdate) {
            await ctx.update({ message: 'Building OBJ mesh', current: t, max: triangleCount });
        }
    }

    const m = MeshBuilder.getMesh(builderState);
    if (!hasNormals) Mesh.computeNormals(m);

    // TODO: check if needed
    ValueCell.updateIfChanged(m.varyingGroup, true);

    return m;
}

/**
 * Build a mesh from indexed OBJ data.
 *
 * OBJ stores positions and normals as indexed arrays (each face-vertex
 * references a position index and an optional normal index). Mesh.create
 * requires a flat vertex array, so we expand the indexed data.
 *
 * When `usePositionGroups` is true (vertex colors present) the group ID for
 * each mesh vertex is set to the position index so the color function can
 * look up per-vertex colors directly.  Otherwise group IDs are material indices
 * (existing behaviour).
 */
async function getMesh(ctx: RuntimeContext, obj: ObjFile, usePositionGroups: boolean, mesh?: Mesh): Promise<Mesh> {
    const { canUseIndexedPath, posGroups, posNormals } = getCanUseIndexedPath(obj, usePositionGroups);

    // Indexed path: map each position directly to a vertex (no duplication).
    // Conditions: (1) no position is shared across different face groups, and
    //             (2) if normals are present, each position maps to exactly one normal index.
    if (canUseIndexedPath) {
        return getIndexedMesh(ctx, obj, { posGroups, posNormals, usePositionGroups }, mesh);
    }

    // Fallback (position uses multiple normals, or material group-conflict detected): expand
    // indexed data into flat per-face-vertex arrays so each face-vertex can carry its own
    // normal/group.
    return getFaceMesh(ctx, obj, usePositionGroups, mesh);
}

function createShape(
    objData: ObjData,
    mesh: Mesh,
    colors: Color[],
    coloringName: ColoringName,
    positionToMaterial: Int32Array | null,
) {
    const { source, transforms } = objData;
    const { materialNames, vertexColors, positionCount } = source;
    const usePositionGroups = vertexColors.length > 0;

    let getColor: (groupId: number) => Color;
    if (coloringName === 'vertex' && usePositionGroups) {
        getColor = (groupId: number) => Color.fromNormalizedArray(vertexColors, groupId * 3);
    } else if (usePositionGroups && positionToMaterial) {
        // material or uniform coloring with position-index groups
        if (coloringName === 'uniform') {
            const c = colors[0] ?? ColorNames.grey;
            getColor = () => c;
        } else {
            getColor = (groupId: number) => colors[Math.min(positionToMaterial[groupId], colors.length - 1)];
        }
    } else {
        getColor = (groupId: number) => colors[Math.min(groupId, colors.length - 1)];
    }

    let getLabel: (groupId: number) => string;
    if (usePositionGroups && positionToMaterial && materialNames.length > 0) {
        getLabel = (groupId: number) => materialNames[positionToMaterial[groupId]] ?? 'OBJ Mesh';
    } else if (!usePositionGroups && materialNames.length > 0) {
        getLabel = (groupId: number) => materialNames[groupId] ?? 'OBJ Mesh';
    } else {
        getLabel = () => 'OBJ Mesh';
    }

    // When using position-index groups the group count equals positionCount;
    // otherwise pass colors.length (existing behaviour — let Shape derive from geometry
    // when colors.length is 0 to avoid a zero group count).
    const groupCount = usePositionGroups ? positionCount : (colors.length > 0 ? colors.length : undefined);

    return Shape.create(
        'obj-mesh', source, mesh,
        (groupId: number) => getColor(groupId),
        () => 1,
        (groupId: number) => getLabel(groupId),
        transforms,
        groupCount
    );
}

function makeShapeGetter() {
    let _objData: ObjData | undefined;
    let _colors: Color[] | undefined;
    let _coloringName: ColoringName | undefined;

    let _shape: Shape<Mesh>;
    let _mesh: Mesh;
    let _positionToMaterial: Int32Array | null = null;

    const getShape = async (ctx: RuntimeContext, objData: ObjData, props: PD.Values<ObjShapeParams>, shape?: Shape<Mesh>) => {
        const newMesh = !_objData || _objData !== objData;
        const coloringName = props.coloring.name;

        const nextColors = getMaterialColors(objData.source.materialNames, props, objData.mtl);
        const newColor = !_colors
            || _coloringName !== coloringName
            || nextColors.length !== _colors.length
            || nextColors.some((c, i) => c !== _colors![i]);

        if (newMesh) {
            _colors = nextColors;
            _coloringName = coloringName;
            const usePositionGroups = objData.source.vertexColors.length > 0;
            _positionToMaterial = usePositionGroups ? buildPositionToMaterial(objData.source) : null;
            _mesh = await getMesh(ctx, objData.source, usePositionGroups, shape && shape.geometry);
            _shape = createShape(objData, _mesh, _colors, _coloringName, _positionToMaterial);
        } else if (newColor) {
            _colors = nextColors;
            _coloringName = coloringName;
            _shape = createShape(objData, _mesh, _colors, _coloringName, _positionToMaterial);
        }

        _objData = objData;

        return _shape;
    };
    return getShape;
}

export function shapeFromObj(source: ObjFile, params?: { transforms?: Mat4[], mtl?: MtlFile }) {
    return Task.create<ShapeProvider<ObjData, Mesh, ObjShapeParams>>('Shape Provider', async _ctx => {
        return {
            label: 'Mesh',
            data: { source, transforms: params?.transforms, mtl: params?.mtl },
            params: createObjShapeParams(source, params?.mtl),
            getShape: makeShapeGetter(),
            geometryUtils: Mesh.Utils,
        };
    });
}
