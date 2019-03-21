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
    if (!nx || !ny || !nz) throw new Error('missing normal properties') // TODO calculate normals when not provided

    for (let i = 0, il = vertex.rowCount; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding vertex ${i}` })

        ChunkedArray.add3(vertices, x.value(i), y.value(i), z.value(i))
        ChunkedArray.add3(normals, nx.value(i), ny.value(i), nz.value(i));
        ChunkedArray.add(groups, groupIds[i])
    }

    for (let i = 0, il = face.rowCount; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding face ${i}` })

        const { entries } = face.value(i)
        ChunkedArray.add3(indices, entries[0], entries[1], entries[2])
    }
    return MeshBuilder.getMesh(builderState);
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

    const red = vertex.getProperty(vp.red)
    const green = vertex.getProperty(vp.green)
    const blue = vertex.getProperty(vp.blue)
    if (!red || !green || !blue) throw new Error('missing color properties')

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
        red: PD.Select('red' as string, [['red', 'red']]),
        green: PD.Select('green' as string, [['green', 'green']]),
        blue: PD.Select('blue' as string, [['blue', 'blue']]),
    }, { isExpanded: true }),
}
export type PlyShapeParams = typeof PlyShapeParams


export function getPlyShapeParams(plyFile: PlyFile) {
    const params = PD.clone(PlyShapeParams)
    const vertex = plyFile.getElement('vertex') as PlyTable
    if (vertex) {
        const options: [string, string][] = [['', '']]
        for (let i = 0, il = vertex.propertyNames.length; i < il; ++i) {
            const name = vertex.propertyNames[i]
            options.push([ name, name ])
        }
        const vp = params.vertexProperties.params;
        (vp.group as PD.Select<string>).options = options;
        (vp.red as PD.Select<string>).options = options;
        (vp.green as PD.Select<string>).options = options;
        (vp.blue as PD.Select<string>).options = options;

        // TODO harcoded as convenience for data provided by MegaMol
        if (vertex.propertyNames.includes('atomid')) {
            vp.group.defaultValue = 'atomid'
            params.vertexProperties.defaultValue.group = 'atomid'
        }
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