/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// TODO

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