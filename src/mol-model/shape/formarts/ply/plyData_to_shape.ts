import {ply_form, PlyFile} from '../../../../mol-io/reader/ply/parse_data/data-model';
import {RuntimeContext, Task} from 'mol-task';
import {Mesh} from '../../../../mol-geo/geometry/mesh/mesh';
import {MeshBuilder} from '../../../../mol-geo/geometry/mesh/mesh-builder';
import { addSphere } from 'mol-geo/geometry/mesh/builder/sphere';
import {Vec3} from '../../../../mol-math/linear-algebra/3d';
import {Shape} from '../../shape';
import {Color} from '../../../../mol-util/color';
import { ShapeProvider } from 'mol-model/shape/provider';

export interface MyData {
    centers: number[],
    colors: Color[],
    labels: string[],
    transforms: number[]
}

function collectData_for_Shape(parsedData: ply_form): MyData {
    // parsedData.data.PLY_File. to access So.format.Ply
    console.log('parsedData', parsedData)
    const { vertices, colors } = parsedData
    const data: MyData = {
        centers: vertices,
        colors: [],
        labels: [],
        transforms: []
    }

    for (let i = 0; i<parsedData.vertexCount; i++) {
        data.colors[i] = Color.fromRgb(colors[i*3+0], colors[i*3+1], colors[i*3+2]);
        data.labels[i] = '';
        data.transforms[i] = 0;
    }
    console.log('data', data);
    return data;
}

async function getSphereMesh(ctx: RuntimeContext, centers: number[], mesh?: Mesh) {
    const builderState = MeshBuilder.createState(centers.length * 128, centers.length * 128 / 2, mesh)
    const v = Vec3.zero()
    builderState.currentGroup = 0
    for (let i = 0, il = centers.length / 3; i < il; ++i) {
        if (i % 10000 === 0 && ctx.shouldUpdate) await ctx.update({ current: i, max: il, message: `adding sphere ${i}` })
        builderState.currentGroup = i
        addSphere(builderState, Vec3.fromArray(v, centers, i * 3), 0.2, 1)
    }
    let a = MeshBuilder.getMesh(builderState);
    // console.log(a);
    return a
}


export async function getShape(ctx: RuntimeContext, parsedData: ply_form, props: {}, shape?: Shape<Mesh>) {
    const data = collectData_for_Shape(parsedData)
    await ctx.update('async creation of shape from  myData')
    const { centers , colors, labels } = data
    const mesh = await getSphereMesh(ctx, centers, shape && shape.geometry)
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