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
import { ColorNames } from 'mol-util/color/tables';

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

    const { coloring, grouping } = props

    const vertex = plyFile.getElement('vertex') as PlyTable
    if (!vertex) throw new Error('missing vertex element')

    const { rowCount } = vertex
    const int = Column.Schema.int

    let red: Column<number>, green: Column<number>, blue: Column<number>
    if (coloring.name === 'vertex') {
        red = vertex.getProperty(coloring.params.red) || Column.ofConst(127, rowCount, int)
        green = vertex.getProperty(coloring.params.green) || Column.ofConst(127, rowCount, int)
        blue = vertex.getProperty(coloring.params.blue) || Column.ofConst(127, rowCount, int)
    } else {
        const [r, g, b] = Color.toRgb(coloring.params.color)
        red = Column.ofConst(r, rowCount, int)
        green = Column.ofConst(g, rowCount, int)
        blue = Column.ofConst(b, rowCount, int)
    }

    const face = plyFile.getElement('face') as PlyList
    if (!face) throw new Error('missing face element')

    const { ids, map } = getGrouping(vertex.rowCount, grouping.name === 'vertex' ? vertex.getProperty(grouping.params.group) : undefined)

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

function createPlyShapeParams(vertex?: PlyTable) {
    const options: [string, string][] = [['', '']]
    const defaultValues = { group: '', red: '', green: '', blue: '' }
    if (vertex) {
        for (let i = 0, il = vertex.propertyNames.length; i < il; ++i) {
            const name = vertex.propertyNames[i]
            options.push([ name, name ])
        }

        // TODO harcoded as convenience for data provided by MegaMol
        if (vertex.propertyNames.includes('atomid')) defaultValues.group = 'atomid'

        if (vertex.propertyNames.includes('red')) defaultValues.red = 'red'
        if (vertex.propertyNames.includes('green')) defaultValues.green = 'green'
        if (vertex.propertyNames.includes('blue')) defaultValues.blue = 'blue'
    }

    return {
        ...Mesh.Params,

        coloring: PD.MappedStatic(defaultValues.red && defaultValues.green && defaultValues.blue ? 'vertex' : 'uniform', {
            vertex: PD.Group({
                red: PD.Select(defaultValues.red, options, { label: 'Red Property' }),
                green: PD.Select(defaultValues.green, options, { label: 'Green Property' }),
                blue: PD.Select(defaultValues.blue, options, { label: 'Blue Property' }),
            }, { isFlat: true }),
            uniform: PD.Group({
                color: PD.Color(ColorNames.grey)
            }, { isFlat: true })
        }),
        grouping: PD.MappedStatic(defaultValues.group ? 'vertex' : 'none', {
            vertex: PD.Group({
                group: PD.Select(defaultValues.group, options, { label: 'Group Property' }),
            }, { isFlat: true }),
            none: PD.Group({ })
        }),
    }
}

export const PlyShapeParams = createPlyShapeParams()
export type PlyShapeParams = typeof PlyShapeParams

export function shapeFromPly(source: PlyFile, params?: {}) {
    return Task.create<ShapeProvider<PlyFile, Mesh, PlyShapeParams>>('Shape Provider', async ctx => {
        return {
            label: 'Mesh',
            data: source,
            params: createPlyShapeParams(source.getElement('vertex') as PlyTable),
            getShape,
            geometryUtils: Mesh.Utils
        }
    })

}