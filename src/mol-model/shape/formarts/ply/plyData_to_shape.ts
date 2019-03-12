import {ply_form, PlyFile} from '../../../../mol-io/reader/ply/parse_data/data-model';
import {RuntimeContext, Task} from 'mol-task';
import {Mesh} from '../../../../mol-geo/geometry/mesh/mesh';
import {MeshBuilder} from '../../../../mol-geo/geometry/mesh/mesh-builder';
import { addTriangle } from 'mol-geo/geometry/mesh/builder/triangle';
import {Shape} from '../../shape';
import {Color} from '../../../../mol-util/color';
import { ShapeProvider } from 'mol-model/shape/provider';

export interface MyData {
    centers: number[],
    normals: number[],
    faces: number[],
    colors: Color[],
    labels: string[],
    transforms: number[]
}

function collectData_for_Shape(parsedData: ply_form): MyData {
    // parsedData.data.PLY_File. to access So.format.Ply
    console.log('parsedData', parsedData)
    const { vertices, colors, faces, normals } = parsedData
    const data: MyData = {
        centers: vertices,
        normals: normals,
        faces: faces,
        colors: [],
        labels: [],
        transforms: []
    }

    for (let i = 0; i<parsedData.faceCount; i++) {
        data.colors[i] = Color.fromRgb(colors[faces[4*i+1]*3+0], colors[faces[4*i+1]*3+1], colors[faces[4*i+1]*3+2]);
        data.labels[i] = parsedData.properties[parsedData.propertyCount * faces[4*i+1] + 10].toString();
            //i.toString();
        // data.transforms[i] = 0;
    }
    console.log('data', data);
    return data;
}

async function getSphereMesh(ctx: RuntimeContext, centers: number[], normals: number[], faces: number[], mesh?: Mesh) {
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
        triangle_indices = [0,1,2];
        //console.log(triangle_vertices)
        addTriangle(builderState, triangle_vertices, triangle_normals, triangle_indices)
    }
    let a = MeshBuilder.getMesh(builderState);
    // a.normalsComputed = false
    // Mesh.computeNormalsImmediate(a)
    console.log(a);
    return a
}



export async function getShape(ctx: RuntimeContext, parsedData: ply_form, props: {}, shape?: Shape<Mesh>) {
    const data = collectData_for_Shape(parsedData)
    await ctx.update('async creation of shape from  myData')
    const { centers, normals, faces, colors, labels } = data
    const mesh = await getSphereMesh(ctx, centers, normals, faces, shape && shape.geometry)
    const groupCount = centers.length / 3
    return shape || Shape.create(
        'test', mesh,
        (groupId: number) => colors[groupId], // color: per group, same for instances
        () => 1, // size: constant
        (groupId: number, instanceId: number) => labels[instanceId * groupCount + groupId] // label: per group and instance
    )
}

export const PlyShapeParams = {
    ...Mesh.Params
}
export type PlyShapeParams = typeof PlyShapeParams

export function shapeFromPly(source: PlyFile, params?: {}) {
    return Task.create<ShapeProvider<ply_form, Mesh, PlyShapeParams>>('Parse Shape Data', async ctx => {
        console.log('source', source)
        return {
            label: 'Mesh',
            data: source.PLY_File,
            getShape,
            geometryUtils: Mesh.Utils
        }
    })
}