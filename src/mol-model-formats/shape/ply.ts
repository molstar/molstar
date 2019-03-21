/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Schäfer, Marco <marco.schaefer@uni-tuebingen.de>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext, Task } from 'mol-task';
import { ShapeProvider } from 'mol-model/shape/provider';
import { Color } from 'mol-util/color';
import { PlyFile, PlyTable, PlyList } from 'mol-io/reader/ply/schema';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Shape } from 'mol-model/shape';
import { ChunkedArray } from 'mol-data/util';
import { arrayMax, fillSerial } from 'mol-util/array';
import { Column } from 'mol-data/db';
import { ParamDefinition as PD } from 'mol-util/param-definition';

// TODO support 'edge' and 'material' elements, see https://www.mathworks.com/help/vision/ug/the-ply-format.html

async function getPlyMesh(ctx: RuntimeContext, vertex: PlyTable, face: PlyList, groupIds: ArrayLike<number>, mesh?: Mesh) {
    const builderState = MeshBuilder.createState(vertex.rowCount, vertex.rowCount / 4, mesh)
    const { vertices, normals, indices, groups } = builderState

    const x = vertex.getProperty('x')
    const y = vertex.getProperty('y')
    const z = vertex.getProperty('z')
    if (!x || !y || !z) throw new Error('missing coordinate properties')

    const nx = vertex.getProperty('nx')
    const ny = vertex.getProperty('ny')
    const nz = vertex.getProperty('nz')

    const hasNormals = !!nx && !!ny && !!nz

    for (let i = 0, il = vertex.rowCount; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding vertex ${i}` })

        ChunkedArray.add3(vertices, x.value(i), y.value(i), z.value(i))
        if (hasNormals) ChunkedArray.add3(normals, nx!.value(i), ny!.value(i), nz!.value(i));
        ChunkedArray.add(groups, groupIds[i])
    }

    for (let i = 0, il = face.rowCount; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding face ${i}` })

        const { entries } = face.value(i)
        ChunkedArray.add3(indices, entries[0], entries[1], entries[2])
    }

    const m = MeshBuilder.getMesh(builderState);
    m.normalsComputed = hasNormals
    await Mesh.computeNormals(m).runInContext(ctx)

    return m
}

function getGrouping(count: number, column?: Column<number>) {
    const ids = column ? column.toArray({ array: Int32Array }) : fillSerial(new Uint32Array(count))
    const maxId = arrayMax(ids) // assumes uint ids
    const map = new Uint32Array(maxId + 1)
    for (let i = 0, il = ids.length; i < il; ++i) map[ids[i]] = i
    return { ids, map }
}

async function getShape(ctx: RuntimeContext, plyFile: PlyFile, props: PD.Values<PlyShapeParams>, shape?: Shape<Mesh>) {
    await ctx.update('async creation of shape from ply file')

    const { vertexProperties: vp } = props

    const vertex = plyFile.getElement('vertex') as PlyTable
    if (!vertex) throw new Error('missing vertex element')

    const red = vertex.getProperty(vp.red) || Column.ofConst(127, vertex.rowCount, Column.Schema.int)
    const green = vertex.getProperty(vp.green) || Column.ofConst(127, vertex.rowCount, Column.Schema.int)
    const blue = vertex.getProperty(vp.blue) || Column.ofConst(127, vertex.rowCount, Column.Schema.int)

    const face = plyFile.getElement('face') as PlyList
    if (!face) throw new Error('missing face element')

    const { ids, map } = getGrouping(vertex.rowCount, vertex.getProperty(vp.group))

    const mesh = await getPlyMesh(ctx, vertex, face, ids, shape && shape.geometry)
    return Shape.create(
        'test', plyFile, mesh,
        (groupId: number) => {
            const idx = map[groupId]
            return Color.fromRgb(red.value(idx), green.value(idx), blue.value(idx))
        },
        () => 1, // size: constant
        (groupId: number) => {
            return ids[groupId].toString()
        }
    )
}

export const PlyShapeParams = {
    ...Mesh.Params,

    vertexProperties: PD.Group({
        group: PD.Select('' as string, [['', '']]),
        red: PD.Select('' as string, [['', '']]),
        green: PD.Select('' as string, [['', '']]),
        blue: PD.Select('' as string, [['', '']]),
    }, { isExpanded: true }),
}
export type PlyShapeParams = typeof PlyShapeParams

function setGroupDefault<T>(group: PD.Group<any>, name: string, defaultValue: T) {
    group.params[name].defaultValue = defaultValue
    group.defaultValue.group = defaultValue
}
function setSelectOptions(select: PD.Select<string>, options: [string, string][]) {
    select.options = options;
}

export function getPlyShapeParams(plyFile: PlyFile) {
    const params = PD.clone(PlyShapeParams)
    const vertex = plyFile.getElement('vertex') as PlyTable
    if (vertex) {
        const options: [string, string][] = [['', '']]
        for (let i = 0, il = vertex.propertyNames.length; i < il; ++i) {
            const name = vertex.propertyNames[i]
            options.push([ name, name ])
        }
        const vp = params.vertexProperties;
        setSelectOptions(vp.params.group as PD.Select<string>, options);
        setSelectOptions(vp.params.red as PD.Select<string>, options);
        setSelectOptions(vp.params.green as PD.Select<string>, options);
        setSelectOptions(vp.params.blue as PD.Select<string>, options);

        // TODO harcoded as convenience for data provided by MegaMol
        if (vertex.propertyNames.includes('atomid')) setGroupDefault(vp, 'group', 'atomid')

        if (vertex.propertyNames.includes('red')) setGroupDefault(vp, 'red', 'red')
        if (vertex.propertyNames.includes('green')) setGroupDefault(vp, 'green', 'green')
        if (vertex.propertyNames.includes('blue')) setGroupDefault(vp, 'blue', 'blue')
    }
    return params
}

export function shapeFromPly(source: PlyFile, params?: {}) {
    return Task.create<ShapeProvider<PlyFile, Mesh, PlyShapeParams>>('Shape Provider', async ctx => {
        return {
            label: 'Mesh',
            data: source,
            params: getPlyShapeParams(source),
            getShape,
            geometryUtils: Mesh.Utils
        }
    })
}