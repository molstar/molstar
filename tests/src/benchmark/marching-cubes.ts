/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Tensor, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { fillSerial } from '@molstar/core/util/array';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { createHistogramPyramid } from '@molstar/graphics/gl/webgl/compute/histogram-pyramid/reduction';
import { calcActiveVoxels } from '@molstar/graphics/gl/webgl/compute/marching-cubes/active-voxels';
import { isComputeIsosurfaceSupported } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import {
  computeChunkedIsosurface,
  destroyChunkedSurface,
  getChunkByteCount,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import { createIsosurfaceBuffers } from '@molstar/graphics/gl/webgl/compute/marching-cubes/isosurface';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';
import { createVolumeTexture2d, getVolumeTexture2dLayout } from '@molstar/graphics/repr/volume/util';
import { Grid, type Volume } from '@molstar/model/model/volume';
import { blobVolume, createVolume } from '../visual/scenes/volume.js';
import { type BenchmarkCase, type PathRun, StageTimer } from './framework.js';

/** Padding of the 2d volume texture, as in the volume isosurface representation */
const Padding = 1;

type VolumeInput = { volume: Volume; isoLevel: number };

/** Gyroid with a period of `period` cells: many surfaces through the whole grid (a dense surface) */
function gyroidVolume(n: number, period: number): Volume {
  const space = Tensor.Space([n, n, n], [2, 1, 0], Float32Array);
  const data = new Float32Array(n * n * n);
  const f = (2 * Math.PI) / period;
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) {
        const x = i * f;
        const y = j * f;
        const z = k * f;
        data[space.dataOffset(i, j, k)] =
          Math.sin(x) * Math.cos(y) + Math.sin(y) * Math.cos(z) + Math.sin(z) * Math.cos(x);
      }
    }
  }
  return createVolume('gyroid', space, data);
}

async function cpuPath(_: GraphicsContext, { volume, isoLevel }: VolumeInput): Promise<PathRun> {
  const timer = new StageTimer();
  const scalarField = volume.grid.cells;
  const ids = fillSerial(new Int32Array(scalarField.data.length));
  const idField = Tensor.create(scalarField.space, Tensor.Data1(ids));
  await timer.stage('setup');
  const mesh = await computeMarchingCubesMesh({ isoLevel, scalarField, idField }).run();
  await timer.stage('marching cubes');
  const bytes =
    mesh.vertexBuffer.ref.value.byteLength +
    mesh.indexBuffer.ref.value.byteLength +
    mesh.normalBuffer.ref.value.byteLength +
    mesh.groupBuffer.ref.value.byteLength;
  return { stages: timer.stages, triangles: mesh.triangleCount, bytes };
}

/** The GPU isosurface of the volume representation (texture upload as there, `gpuDataType` 'byte') */
async function webglPath(graphics: GraphicsContext, { volume, isoLevel }: VolumeInput): Promise<PathRun | string> {
  if (!isWebGLContext(graphics)) return 'needs a WebGL context';
  if (!graphics.capabilities.multipleRenderTargets) return 'needs multiple render targets';
  const timer = new StageTimer(graphics);

  const gridDimension = Vec3.clone(volume.grid.cells.space.dimensions as Vec3);
  const { width, height, powerOfTwoSize: texDim } = getVolumeTexture2dLayout(gridDimension, Padding);
  if (texDim > graphics.resources.limits.maxTextureSize / 2) return `needs a ${texDim}² texture`;
  const gridTexDim = Vec3.create(width, height, 0);
  const gridDataDim = Vec3.clone(gridDimension);
  const gridTexScale = Vec2.create(width / texDim, height / texDim);
  const texture = graphics.resources.texture('image-uint8', 'r', 'ubyte', 'linear');
  texture.define(texDim, texDim);
  texture.load(createVolumeTexture2d(volume, 'data', Padding, 'byte'), true);
  gridDimension[0] += Padding;
  gridDimension[1] += Padding;
  await timer.stage('upload');

  const { min, max } = volume.grid.stats;
  const level = (isoLevel - min) / (max - min);
  const activeVoxels = calcActiveVoxels(graphics, texture, 'red', gridDimension, gridTexDim, level, gridTexScale);
  await timer.stage('active voxels');
  const pyramid = createHistogramPyramid(graphics, activeVoxels, gridTexScale, gridTexDim);
  await timer.stage('histogram pyramid');
  const transform = Grid.getGridToCartesianTransform(volume.grid);
  const axisOrder = volume.grid.cells.space.axisOrderSlowToFast as Vec3;
  const buffers = createIsosurfaceBuffers(
    graphics,
    activeVoxels,
    texture,
    'red',
    pyramid,
    gridDimension,
    gridTexDim,
    gridDataDim,
    transform,
    level,
    isoLevel < 0,
    false,
    axisOrder,
    true,
  );
  await timer.stage('isosurface');

  const bytes =
    buffers.vertexTexture.getByteCount() + buffers.groupTexture.getByteCount() + buffers.normalTexture.getByteCount();
  const triangles = buffers.vertexCount / 3;
  texture.destroy();
  buffers.vertexTexture.destroy();
  buffers.groupTexture.destroy();
  buffers.normalTexture.destroy();
  return { stages: timer.stages, triangles, bytes };
}

/** Marching cubes with compute kernels on the raw values in a storage buffer, into the vertex arena; until the count is read back */
async function webgpuPath(graphics: GraphicsContext, { volume, isoLevel }: VolumeInput): Promise<PathRun | string> {
  if (!isComputeIsosurfaceSupported(graphics)) return 'needs compute shaders';
  const timer = new StageTimer(graphics);
  const { space, data } = volume.grid.cells;
  const values = data instanceof Float32Array ? data : Float32Array.from(data as ArrayLike<number>);
  const buffer = graphics.compute.createBuffer(values.byteLength, 'benchmark volume');
  buffer.write(values);
  const [sx, sy, sz] = [0, 1, 2].map((axis) =>
    space.dataOffset(axis === 0 ? 1 : 0, axis === 1 ? 1 : 0, axis === 2 ? 1 : 0),
  );
  await timer.stage('upload');
  const transform = Grid.getGridToCartesianTransform(volume.grid);
  const { compute } = graphics;
  const result = await computeChunkedIsosurface(
    compute,
    { kind: 'buffer', volume: buffer, strides: Vec3.create(sx, sy, sz) },
    Vec3.clone(space.dimensions as Vec3),
    isoLevel,
    transform,
    isoLevel < 0,
  );
  const { count, pages } = await result.count;
  await timer.stage('marching cubes');
  // the bytes of the pages taken
  const bytes = getChunkByteCount(pages);
  buffer.destroy();
  destroyChunkedSurface(compute, result.surface);
  return { stages: timer.stages, triangles: count / 3, bytes };
}

function volumeCase(id: string, description: string, setup: () => VolumeInput): BenchmarkCase<VolumeInput> {
  return { id: `mc/volume/${id}`, description, setup, paths: { cpu: cpuPath, webgl: webglPath, webgpu: webgpuPath } };
}

export const MarchingCubesCases: BenchmarkCase<any>[] = [64, 128, 192, 256].flatMap((n) => [
  volumeCase(`blobs-${n}`, `Three Gaussian blobs on ${n}³ cells, iso 0.5: a sparse surface.`, () => ({
    volume: blobVolume(n),
    isoLevel: 0.5,
  })),
  volumeCase(`gyroid-${n}`, `Gyroid with a period of 16 cells on ${n}³ cells, iso 0: a dense surface.`, () => ({
    volume: gyroidVolume(n, 16),
    isoLevel: 0,
  })),
]);
