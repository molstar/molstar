/**
 * Copyright (c) 2019 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Schäfer, Marco <marco.schaefer@uni-tuebingen.de>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { RuntimeContext, Task } from 'mol-task';
import { addTriangle } from 'mol-geo/geometry/mesh/builder/triangle';
import { ShapeProvider } from 'mol-model/shape/provider';
import { Color } from 'mol-util/color';
import { PlyData, PlyFile } from 'mol-io/reader/ply/schema';
import { MeshBuilder } from 'mol-geo/geometry/mesh/mesh-builder';
import { Mesh } from 'mol-geo/geometry/mesh/mesh';
import { Shape } from 'mol-model/shape';

async function getPlyMesh(ctx: RuntimeContext, centers: number[], normals: number[], faces: number[], mesh?: Mesh) {
    const builderState = MeshBuilder.createState(faces.length, faces.length, mesh)
    builderState.currentGroup = 0
    for (let i = 0, il = faces.length/4; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding triangle ${i}` })
        builderState.currentGroup = i

        let triangle_vertices: number[];
        let triangle_normals: number[];
        let triangle_indices: number[];
        triangle_vertices = [centers[faces[4*i+1]*3], centers[faces[4*i+1]*3+1], centers[faces[4*i+1]*3+2],
                             centers[faces[4*i+2]*3], centers[faces[4*i+2]*3+1], centers[faces[4*i+2]*3+2],
                             centers[faces[4*i+3]*3], centers[faces[4*i+3]*3+1], centers[faces[4*i+3]*3+2]];
        triangle_normals = [ normals[faces[4*i+1]*3], normals[faces[4*i+1]*3+1], normals[faces[4*i+1]*3+2],
                             normals[faces[4*i+2]*3], normals[faces[4*i+2]*3+1], normals[faces[4*i+2]*3+2],
                             normals[faces[4*i+3]*3], normals[faces[4*i+3]*3+1], normals[faces[4*i+3]*3+2]];
        triangle_indices = [0, 1, 2];
        // console.log(triangle_vertices)
        addTriangle(builderState, triangle_vertices, triangle_normals, triangle_indices)
    }
    return MeshBuilder.getMesh(builderState);
}

async function getShape(ctx: RuntimeContext, parsedData: PlyData, props: {}, shape?: Shape<Mesh>) {
    await ctx.update('async creation of shape from  myData')
    const { vertices, normals, faces, colors, properties } = parsedData
    const mesh = await getPlyMesh(ctx, vertices, normals, faces, shape && shape.geometry)
    return shape || Shape.create(
        'ply-mesh', mesh,
        (groupId: number) => {
            return Color.fromRgb(
                colors[faces[4 * groupId + 1] * 3 + 0],
                colors[faces[4 * groupId + 1] * 3 + 1],
                colors[faces[4 * groupId + 1] * 3 + 2]
            )
        },
        () => 1, // size: constant
        (groupId: number) => {
            return properties[parsedData.propertyCount * faces[4 * groupId + 1] + 10].toString()
        }
    )
}

export const PlyShapeParams = {
    ...Mesh.Params
}
export type PlyShapeParams = typeof PlyShapeParams

export function shapeFromPly(source: PlyFile, params?: {}) {
    return Task.create<ShapeProvider<PlyData, Mesh, PlyShapeParams>>('Parse Shape Data', async ctx => {
        console.log('source', source)
        return {
            label: 'Mesh',
            data: source.PLY_File,
            getShape,
            geometryUtils: Mesh.Utils
        }
    })

}