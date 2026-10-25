/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Vec3 } from '@molstar/core/math/linear-algebra';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { isComputeIsosurfaceSupported } from '@molstar/graphics/gl/compute/marching-cubes/compute';
import {
  computeChunkedIsosurface,
  destroyChunkedSurface,
  getChunkByteCount,
} from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import { extractIsosurface } from '@molstar/graphics/gl/webgl/compute/marching-cubes/isosurface';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';
import {
  computeStructureGaussianDensity,
  computeStructureGaussianDensityTexture,
  computeStructureGaussianDensityTexture2d,
  GaussianDensityParams,
} from '@molstar/graphics/repr/structure/visual/util/gaussian';
import type { SizeTheme } from '@molstar/graphics/theme/size';
import type { Structure } from '@molstar/model/model/structure';
import { type BenchmarkCase, type PathRun, StageTimer } from './framework.js';
import { exampleStructure, getPhysicalSizeTheme, globuleStructure } from './structures.js';

type GaussianInput = {
  structure: Structure;
  sizeTheme: SizeTheme<any>;
  props: PD.Values<typeof GaussianDensityParams>;
};

/** The structure Gaussian surface mesh visual: CPU density, then CPU marching cubes */
async function cpuPath(_: GraphicsContext, { structure, sizeTheme, props }: GaussianInput): Promise<PathRun> {
  const timer = new StageTimer();
  const { field, idField, radiusFactor } = await computeStructureGaussianDensity(structure, sizeTheme, props).run();
  await timer.stage('density');
  const isoLevel = Math.exp(-props.smoothness) / radiusFactor;
  const mesh = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }).run();
  await timer.stage('marching cubes');
  const bytes =
    mesh.vertexBuffer.ref.value.byteLength +
    mesh.indexBuffer.ref.value.byteLength +
    mesh.normalBuffer.ref.value.byteLength +
    mesh.groupBuffer.ref.value.byteLength;
  return { stages: timer.stages, triangles: mesh.triangleCount, bytes };
}

/** The density texture, reused across runs as the representation reuses its named texture */
const densityTextures = new WeakMap<GraphicsContext, Texture>();
function getDensityTexture(graphics: GraphicsContext) {
  let texture = densityTextures.get(graphics);
  if (!texture) {
    const { resources, capabilities } = graphics;
    texture = capabilities.halfFloatRenderTarget
      ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
      : capabilities.floatBlend
        ? resources.texture('image-float32', 'rgba', 'float', 'linear')
        : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    densityTextures.set(graphics, texture);
  }
  return texture;
}

/** The structure Gaussian surface texture mesh visual: GPU density (2d atlas), then GPU marching cubes */
async function webglPath(
  graphics: GraphicsContext,
  { structure, sizeTheme, props }: GaussianInput,
): Promise<PathRun | string> {
  if (!isWebGLContext(graphics)) return 'needs a WebGL context';
  if (!graphics.capabilities.multipleRenderTargets) return 'needs multiple render targets';
  const timer = new StageTimer(graphics);
  const density = computeStructureGaussianDensityTexture2d(
    structure,
    sizeTheme,
    true,
    props,
    graphics,
    getDensityTexture(graphics),
  );
  await timer.stage('density');
  const isoLevel = Math.exp(-props.smoothness) / density.radiusFactor;
  const buffers = extractIsosurface(
    graphics,
    density.texture,
    'alpha',
    density.gridDim,
    density.gridTexDim,
    density.gridDataDim,
    density.gridTexScale,
    density.transform,
    isoLevel,
    false,
    true,
    Vec3.create(0, 1, 2),
    true,
  );
  await timer.stage('marching cubes');
  const bytes =
    buffers.vertexTexture.getByteCount() + buffers.groupTexture.getByteCount() + buffers.normalTexture.getByteCount();
  const triangles = buffers.vertexCount / 3;
  buffers.vertexTexture.destroy();
  buffers.groupTexture.destroy();
  buffers.normalTexture.destroy();
  return { stages: timer.stages, triangles, bytes };
}

/** The 3d density texture of the WebGPU path, reused across runs */
const densityTextures3d = new WeakMap<GraphicsContext, Texture>();

/** GPU density (3d texture), then marching cubes with compute kernels into the vertex arena, until the count is read back */
async function webgpuPath(
  graphics: GraphicsContext,
  { structure, sizeTheme, props }: GaussianInput,
): Promise<PathRun | string> {
  if (!isComputeIsosurfaceSupported(graphics)) return 'needs compute shaders';
  const timer = new StageTimer(graphics);
  let texture = densityTextures3d.get(graphics);
  if (!texture) {
    texture = graphics.resources.texture('volume-float16', 'rgba', 'fp16', 'linear');
    densityTextures3d.set(graphics, texture);
  }
  const density = computeStructureGaussianDensityTexture(structure, sizeTheme, props, graphics, texture);
  await timer.stage('density');
  const isoLevel = Math.exp(-props.smoothness) / density.radiusFactor;
  const { compute } = graphics;
  const result = await computeChunkedIsosurface(
    compute,
    { kind: 'texture', texture: density.texture, channel: 'alpha', packedGroup: true },
    density.gridDim,
    isoLevel,
    density.transform,
    false,
  );
  const { count, pages } = await result.count;
  await timer.stage('marching cubes');
  // the bytes of the pages taken
  const bytes = getChunkByteCount(pages);
  destroyChunkedSurface(compute, result.surface);
  return { stages: timer.stages, triangles: count / 3, bytes };
}

function gaussianCase(
  name: string,
  description: string,
  structure: () => Promise<Structure>,
  resolution: number,
  group = 'mc/gaussian',
): BenchmarkCase<GaussianInput> {
  return {
    id: `${group}/${name}-res${resolution}`,
    description: `Gaussian surface of ${description} at a resolution of ${resolution} Å.`,
    setup: async () => {
      const s = await structure();
      const props = { ...PD.getDefaultValues(GaussianDensityParams), resolution };
      return { structure: s, sizeTheme: getPhysicalSizeTheme(s), props };
    },
    paths: { cpu: cpuPath, webgl: webglPath, webgpu: webgpuPath },
  };
}

function globuleCase(atoms: number, resolution: number, group?: string) {
  const label = atoms >= 1000 ? `${atoms / 1000}k` : `${atoms}`;
  const description = `a protein-like ball of ${atoms} atoms`;
  return gaussianCase(`globule-${label}`, description, () => globuleStructure(atoms), resolution, group);
}

function entryCase(entry: string, file: string, description: string, resolution: number) {
  return gaussianCase(
    entry,
    `PDB entry ${entry} (${description}), all units`,
    () => exampleStructure(file),
    resolution,
  );
}

export const GaussianSurfaceCases: BenchmarkCase<any>[] = [
  ...[1000, 10000, 100000].flatMap((atoms) => [globuleCase(atoms, 1), globuleCase(atoms, 0.5)]),
  ...[1, 0.5].flatMap((resolution) => [
    entryCase('1cbs', '1cbs_updated.cif', 'about 1.2 k atoms, a protein with ligand and waters', resolution),
    entryCase('7qpd', '7qpd.fw2.cif', 'about 11.5 k atoms in 7 chains', resolution),
  ]),
];

/**
 * Where the CPU and the GPU paths cross (`getGpuFallbackReason`): small structures to large ones, at resolutions on both
 * sides of the CPU cutoff; not in the default runs (`--filter 'mc/gaussian-cutoff/*'`)
 */
export const GaussianSurfaceCutoffCases: BenchmarkCase<any>[] = [100, 300, 1000, 3000, 10000, 30000].flatMap((atoms) =>
  [0.5, 1, 1.5, 2, 3].map((resolution) => globuleCase(atoms, resolution, 'mc/gaussian-cutoff')),
);
