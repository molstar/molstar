/**
 * Copyright (c) 2018-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// Marching cubes into meshes and lines: the computation is core's (`@molstar/core/math/geometry/marching-cubes`),
// the geometry builders graphics' (`builder.ts`)

import { Task, RuntimeContext } from '@molstar/core/task';
import { Tensor } from '@molstar/core/math/linear-algebra';
import { anchorMap } from '@molstar/core/util/map-anchor';
import {
  computeMarchingCubes,
  createMarchingCubesComputation,
  getMarchingCubesInputParams,
  getMarchingCubesVertexChunkSize,
  type MarchingCubesParams,
} from '@molstar/core/math/geometry/marching-cubes/algorithm';
import type { Mesh } from '../../geometry/mesh/mesh.js';
import type { Lines } from '../../geometry/lines/lines.js';
import { MarchingCubesMeshBuilder, MarchingCubesLinesBuilder } from './builder.js';

export type { MarchingCubesParams };

export function computeMarchingCubesMesh(params: MarchingCubesParams, mesh?: Mesh) {
  return Task.create('Marching Cubes Mesh', async (ctx) => {
    const inputParams = getMarchingCubesInputParams(params);
    const builder = MarchingCubesMeshBuilder(getMarchingCubesVertexChunkSize(inputParams), mesh);
    return computeMarchingCubes(ctx, builder, inputParams);
  });
}

export function computeMarchingCubesLines(params: MarchingCubesParams, lines?: Lines) {
  return Task.create('Marching Cubes Lines', async (ctx) => {
    const inputParams = getMarchingCubesInputParams(params);
    const builder = MarchingCubesLinesBuilder(getMarchingCubesVertexChunkSize(inputParams), lines);
    return computeMarchingCubes(ctx, builder, inputParams);
  });
}

// keeps the maps of a computation, its state and the mesh builder, which the optimized code of `slice`, `processCell`,
// `interpolate` & co. embeds (`anchorMap`); made, not run, so that it adds no feedback (perf-audit-3 08/P15)
{
  const space = Tensor.Space([2, 2, 2], [0, 1, 2], Float32Array);
  const params = getMarchingCubesInputParams({ isoLevel: 0.5, scalarField: Tensor.create(space, space.create()) });
  anchorMap(createMarchingCubesComputation(RuntimeContext.Synchronous, MarchingCubesMeshBuilder(1), params));
}
