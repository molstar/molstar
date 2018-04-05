/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { Task } from 'mol-task'

export interface Surface {
    vertexCount: number,
    triangleCount: number,
    vertexBuffer: Float32Array,
    indexBuffer: Uint32Array,
    normalBuffer?: Float32Array,
    normalsComputed: boolean,

    vertexAnnotation?: ArrayLike<number>[]
    //boundingSphere?: { center: Geometry.LinearAlgebra.Vector3, radius: number };
}

export namespace Surface {
    export function computeNormalsImmediate(surface: Surface) {
        if (surface.normalsComputed) return;

        const normals = surface.normalBuffer && surface.normalBuffer.length >= surface.vertexCount * 3
            ? surface.normalBuffer : new Float32Array(surface.vertexBuffer.length);

        const v = surface.vertexBuffer, triangles = surface.indexBuffer;
        for (let i = 0, ii = 3 * surface.triangleCount; i < ii; i += 3) {
            const a = 3 * triangles[i],
                b = 3 * triangles[i + 1],
                c = 3 * triangles[i + 2];

            const nx = v[a + 2] * (v[b + 1] - v[c + 1]) + v[b + 2] * v[c + 1] - v[b + 1] * v[c + 2] + v[a + 1] * (-v[b + 2] + v[c + 2]),
                ny = -(v[b + 2] * v[c]) + v[a + 2] * (-v[b] + v[c]) + v[a] * (v[b + 2] - v[c + 2]) + v[b] * v[c + 2],
                nz = v[a + 1] * (v[b] - v[c]) + v[b + 1] * v[c] - v[b] * v[c + 1] + v[a] * (-v[b + 1] + v[b + 1]);

            normals[a] += nx; normals[a + 1] += ny; normals[a + 2] += nz;
            normals[b] += nx; normals[b + 1] += ny; normals[b + 2] += nz;
            normals[c] += nx; normals[c + 1] += ny; normals[c + 2] += nz;
        }

        for (let i = 0, ii = 3 * surface.vertexCount; i < ii; i += 3) {
            const nx = normals[i];
            const ny = normals[i + 1];
            const nz = normals[i + 2];
            const f = 1.0 / Math.sqrt(nx * nx + ny * ny + nz * nz);
            normals[i] *= f; normals[i + 1] *= f; normals[i + 2] *= f;
        }
        surface.normalBuffer = normals;
        surface.normalsComputed = true;
    }

    export function computeNormals(surface: Surface): Task<Surface> {
        return Task.create<Surface>('Surface (Compute Normals)', async ctx => {
            if (surface.normalsComputed) return surface;

            await ctx.update('Computing normals...');
            computeNormalsImmediate(surface);
            return surface;
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

//     export function computeBoundingSphere(surface: Surface): Computation<Surface> {
//         return computation<Surface>(async ctx => {
//             if (surface.boundingSphere) {
//                 return surface;
//             }
//             await ctx.updateProgress('Computing bounding sphere...');

//             const vertices = surface.vertices;
//             let x = 0, y = 0, z = 0;
//             for (let i = 0, _c = surface.vertices.length; i < _c; i += 3) {
//                 x += vertices[i];
//                 y += vertices[i + 1];
//                 z += vertices[i + 2];
//             }
//             x /= surface.vertexCount;
//             y /= surface.vertexCount;
//             z /= surface.vertexCount;
//             let r = 0;
//             for (let i = 0, _c = vertices.length; i < _c; i += 3) {
//                 const dx = x - vertices[i];
//                 const dy = y - vertices[i + 1];
//                 const dz = z - vertices[i + 2];
//                 r = Math.max(r, dx * dx + dy * dy + dz * dz);
//             }
//             surface.boundingSphere = {
//                 center: LinearAlgebra.Vector3.fromValues(x, y, z),
//                 radius: Math.sqrt(r)
//             }
//             return surface;
//         });
//     }

//     export function transformImmediate(surface: Surface, t: number[]) {
//         const p = LinearAlgebra.Vector3.zero();
//         const m = LinearAlgebra.Vector3.transformMat4;
//         const vertices = surface.vertices;
//         for (let i = 0, _c = surface.vertices.length; i < _c; i += 3) {
//             p[0] = vertices[i];
//             p[1] = vertices[i + 1];
//             p[2] = vertices[i + 2];
//             m(p, p, t);
//             vertices[i] = p[0];
//             vertices[i + 1] = p[1];
//             vertices[i + 2] = p[2];
//         }
//         surface.normals = void 0;
//         surface.boundingSphere = void 0;
//     }

//     export function transform(surface: Surface, t: number[]): Computation<Surface> {
//         return computation<Surface>(async ctx => {
//             ctx.updateProgress('Updating surface...');
//             transformImmediate(surface, t);
//             return surface;
//         });
//     }
// }