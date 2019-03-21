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

async function getPlyMesh(ctx: RuntimeContext, vertex: PlyTable, face: PlyList, mesh?: Mesh) {
    const builderState = MeshBuilder.createState(face.rowCount, face.rowCount, mesh)
    const { vertices, normals, indices, groups } = builderState

    const x = vertex.getProperty('x')
    const y = vertex.getProperty('y')
    const z = vertex.getProperty('z')
    if (!x || !y || !z) throw new Error('missing coordinate properties')

    const nx = vertex.getProperty('nx')
    const ny = vertex.getProperty('ny')
    const nz = vertex.getProperty('nz')
    if (!nx || !ny || !nz) throw new Error('missing normal properties')

    const atomid = vertex.getProperty('atomid')
    if (!atomid) throw new Error('missing atomid property')

    for (let i = 0, il = vertex.rowCount; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding vertex ${i}` })

        ChunkedArray.add3(vertices, x.value(i), y.value(i), z.value(i))
        ChunkedArray.add3(normals, nx.value(i), ny.value(i), nz.value(i));
        ChunkedArray.add(groups, atomid.value(i))
    }

    for (let i = 0, il = face.rowCount; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding face ${i}` })

        const { entries } = face.value(i)
        ChunkedArray.add3(indices, entries[0], entries[1], entries[2])
    }
    return MeshBuilder.getMesh(builderState);
}

async function getShape(ctx: RuntimeContext, plyFile: PlyFile, props: {}, shape?: Shape<Mesh>) {
    await ctx.update('async creation of shape from  myData')

    const vertex = plyFile.getElement('vertex') as PlyTable
    if (!vertex) throw new Error('missing vertex element')

    const atomid = vertex.getProperty('atomid')
    if (!atomid) throw new Error('missing atomid property')

    const red = vertex.getProperty('red')
    const green = vertex.getProperty('green')
    const blue = vertex.getProperty('blue')
    if (!red || !green || !blue) throw new Error('missing color properties')

    const face = plyFile.getElement('face') as PlyList
    if (!face) throw new Error('missing face element')

    const mesh = await getPlyMesh(ctx, vertex, face, shape && shape.geometry)
    return shape || Shape.create(

        'test', plyFile, mesh,

        (groupId: number) => {
            return Color.fromRgb(red.value(groupId), green.value(groupId), blue.value(groupId))
        },
        () => 1, // size: constant
        (groupId: number) => {
            return atomid.value(groupId).toString()
        }
    )
}

export const PlyShapeParams = {
    ...Mesh.Params
}
export type PlyShapeParams = typeof PlyShapeParams

export function shapeFromPly(source: PlyFile, params?: {}) {
    return Task.create<ShapeProvider<PlyFile, Mesh, PlyShapeParams>>('Parse Shape Data', async ctx => {
        console.log('source', source)
        return {
            label: 'Mesh',
            data: source,
            getShape,
            geometryUtils: Mesh.Utils
        }
    })

}