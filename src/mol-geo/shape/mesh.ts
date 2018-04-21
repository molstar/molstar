/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from 'mol-task'
import { ValueCell } from 'mol-util'
import { Vec3, Mat4 } from 'mol-math/linear-algebra'
import Sphere from 'mol-math/geometry/sphere'
import { transformPositionArray } from '../util';

export interface Mesh {
    /** Number of vertices in the mesh */
    vertexCount: number,
    /** Number of triangles in the mesh */
    triangleCount: number,

    /** Vertex buffer as array of xyz values wrapped in a value cell */
    vertexBuffer: ValueCell<Float32Array>,
    /** Index buffer as array of vertex index triplets wrapped in a value cell */
    indexBuffer: ValueCell<Uint32Array>,
    /** Normal buffer as array of xyz values for each vertex wrapped in a value cell */
    normalBuffer: ValueCell<Float32Array | undefined>,
    /** Id buffer as array of ids for each vertex wrapped in a value cell */
    idBuffer: ValueCell<Float32Array | undefined>,

    /** Flag indicating if normals are computed for the current set of vertices */
    normalsComputed: boolean,

    /** Bounding sphere of the mesh */
    boundingSphere?: Sphere
}

export namespace Mesh {
    export function computeNormalsImmediate(surface: Mesh) {
        if (surface.normalsComputed) return;

        const normals = surface.normalBuffer.ref.value && surface.normalBuffer.ref.value.length >= surface.vertexCount * 3
            ? surface.normalBuffer.ref.value : new Float32Array(surface.vertexBuffer.ref.value.length);

        const v = surface.vertexBuffer.ref.value, triangles = surface.indexBuffer.ref.value;

        const x = Vec3.zero(), y = Vec3.zero(), z = Vec3.zero(), d1 = Vec3.zero(), d2 = Vec3.zero(), n = Vec3.zero();
        for (let i = 0, ii = 3 * surface.triangleCount; i < ii; i += 3) {
            const a = 3 * triangles[i], b = 3 * triangles[i + 1], c = 3 * triangles[i + 2];

            Vec3.fromArray(x, v, a);
            Vec3.fromArray(y, v, b);
            Vec3.fromArray(z, v, c);
            Vec3.sub(d1, z, y);
            Vec3.sub(d2, y, x);
            Vec3.cross(n, d1, d2);

            normals[a] += n[0]; normals[a + 1] += n[1]; normals[a + 2] += n[2];
            normals[b] += n[0]; normals[b + 1] += n[1]; normals[b + 2] += n[2];
            normals[c] += n[0]; normals[c + 1] += n[1]; normals[c + 2] += n[2];
        }

        for (let i = 0, ii = 3 * surface.vertexCount; i < ii; i += 3) {
            const nx = normals[i];
            const ny = normals[i + 1];
            const nz = normals[i + 2];
            const f = 1.0 / Math.sqrt(nx * nx + ny * ny + nz * nz);
            normals[i] *= f; normals[i + 1] *= f; normals[i + 2] *= f;

            // console.log([normals[i], normals[i + 1], normals[i + 2]], [v[i], v[i + 1], v[i + 2]])
        }
        surface.normalBuffer = ValueCell.update(surface.normalBuffer, normals);
        surface.normalsComputed = true;
    }

    export function computeNormals(surface: Mesh): Task<Mesh> {
        return Task.create<Mesh>('Surface (Compute Normals)', async ctx => {
            if (surface.normalsComputed) return surface;

            await ctx.update('Computing normals...');
            computeNormalsImmediate(surface);
            return surface;
        });
    }

    export function transformImmediate(mesh: Mesh, t: Mat4) {
        transformRangeImmediate(mesh, t, 0, mesh.vertexCount)
    }

    export function transformRangeImmediate(mesh: Mesh, t: Mat4, offset: number, count: number) {
        transformPositionArray(t, mesh.vertexBuffer.ref.value, offset, count)
        // transformDirectionArray(n, mesh.normalBuffer.ref.value, offset, count)  // TODO
        mesh.normalsComputed = false;
        // mesh.boundingSphere = void 0;
    }

    export function computeBoundingSphere(mesh: Mesh): Task<Mesh> {
        return Task.create<Mesh>('Mesh (Compute Bounding Sphere)', async ctx => {
            if (mesh.boundingSphere) {
                return mesh;
            }
            await ctx.update('Computing bounding sphere...');

            const vertices = mesh.vertexBuffer.ref.value;
            let x = 0, y = 0, z = 0;
            for (let i = 0, _c = vertices.length; i < _c; i += 3) {
                x += vertices[i];
                y += vertices[i + 1];
                z += vertices[i + 2];
            }
            x /= mesh.vertexCount;
            y /= mesh.vertexCount;
            z /= mesh.vertexCount;
            let r = 0;
            for (let i = 0, _c = vertices.length; i < _c; i += 3) {
                const dx = x - vertices[i];
                const dy = y - vertices[i + 1];
                const dz = z - vertices[i + 2];
                r = Math.max(r, dx * dx + dy * dy + dz * dz);
            }
            mesh.boundingSphere = {
                center: Vec3.create(x, y, z),
                radius: Math.sqrt(r)
            }
            return mesh;
        });
    }
}

//     function addVertex(src: Float32Array, i: number, dst: Float32Array, j: number) {
//         dst[3 * j] += src[3 * i];
//         dst[3 * j + 1] += src[3 * i + 1];
//         dst[3 * j + 2] += src[3 * i + 2];
//     }

//     function laplacianSmoothIter(surface: Surface, vertexCounts: Int32Array, vs: Float32Array, vertexWeight: number) {
//         const triCount = surface.triangleIndices.length,
//             src = surface.vertices;

//         const triangleIndices = surface.triangleIndices;

//         for (let i = 0; i < triCount; i += 3) {
//             const a = triangleIndices[i],
//                 b = triangleIndices[i + 1],
//                 c = triangleIndices[i + 2];

//             addVertex(src, b, vs, a);
//             addVertex(src, c, vs, a);

//             addVertex(src, a, vs, b);
//             addVertex(src, c, vs, b);

//             addVertex(src, a, vs, c);
//             addVertex(src, b, vs, c);
//         }

//         const vw = 2 * vertexWeight;
//         for (let i = 0, _b = surface.vertexCount; i < _b; i++) {
//             const n = vertexCounts[i] + vw;
//             vs[3 * i] = (vs[3 * i] + vw * src[3 * i]) / n;
//             vs[3 * i + 1] = (vs[3 * i + 1] + vw * src[3 * i + 1]) / n;
//             vs[3 * i + 2] = (vs[3 * i + 2] + vw * src[3 * i + 2]) / n;
//         }
//     }

//     async function laplacianSmoothComputation(ctx: Computation.Context, surface: Surface, iterCount: number, vertexWeight: number) {
//         await ctx.updateProgress('Smoothing surface...', true);

//         const vertexCounts = new Int32Array(surface.vertexCount),
//             triCount = surface.triangleIndices.length;

//         const tris = surface.triangleIndices;
//         for (let i = 0; i < triCount; i++) {
//             // in a triangle 2 edges touch each vertex, hence the constant.
//             vertexCounts[tris[i]] += 2;
//         }

//         let vs = new Float32Array(surface.vertices.length);
//         let started = Utils.PerformanceMonitor.currentTime();
//         await ctx.updateProgress('Smoothing surface...', true);
//         for (let i = 0; i < iterCount; i++) {
//             if (i > 0) {
//                 for (let j = 0, _b = vs.length; j < _b; j++) vs[j] = 0;
//             }
//             surface.normals = void 0;
//             laplacianSmoothIter(surface, vertexCounts, vs, vertexWeight);
//             const t = surface.vertices;
//             surface.vertices = <any>vs;
//             vs = <any>t;

//             const time = Utils.PerformanceMonitor.currentTime();
//             if (time - started > Computation.UpdateProgressDelta) {
//                 started = time;
//                 await ctx.updateProgress('Smoothing surface...', true, i + 1, iterCount);
//             }
//         }
//         return surface;
//     }

//     /*
//      * Smooths the vertices by averaging the neighborhood.
//      *
//      * Resets normals. Might replace vertex array.
//      */
//     export function laplacianSmooth(surface: Surface, iterCount: number = 1, vertexWeight: number = 1): Computation<Surface> {

//         if (iterCount < 1) iterCount = 0;
//         if (iterCount === 0) return Computation.resolve(surface);

//         return computation(async ctx => await laplacianSmoothComputation(ctx, surface, iterCount, (1.1 * vertexWeight) / 1.1));
//     }