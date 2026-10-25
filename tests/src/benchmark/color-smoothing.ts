/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4 } from '@molstar/core/math/linear-algebra';
import { packIntToRGBArray } from '@molstar/core/util/number-packing';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { calcMeshColorSmoothing } from '@molstar/graphics/geo/geometry/mesh/color-smoothing';
import { Mesh } from '@molstar/graphics/geo/geometry/mesh/mesh';
import { calcTextureMeshColorSmoothing } from '@molstar/graphics/geo/geometry/texture-mesh/color-smoothing';
import { computeMarchingCubesMesh } from '@molstar/graphics/geo/util/marching-cubes/algorithm';
import { isComputeColorSmoothingSupported } from '@molstar/graphics/gl/compute/color-smoothing';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { createDataArray, type DataArray, getDataTextureImage } from '@molstar/graphics/gl/renderable/util';
import { createNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { isWebGLContext } from '@molstar/graphics/gl/webgl/context';
import {
  computeStructureGaussianDensity,
  computeStructureGaussianDensityTexture,
  GaussianDensityParams,
} from '@molstar/graphics/repr/structure/visual/util/gaussian';
import type { Structure } from '@molstar/model/model/structure';
import { type ChunkedSurface, computeChunkedIsosurface } from '@molstar/graphics/gl/compute/marching-cubes/chunks';
import { type BenchmarkCase, type PathRun, StageTimer } from './framework.js';
import { getPhysicalSizeTheme, globuleStructure } from './structures.js';

/** Resolution of the surfaces, Å */
const SurfaceResolution = 0.5;
/** Every 3rd vertex is sampled, the default `sampleStride` */
const SampleStride = 3;

type SmoothingInput = {
  structure: Structure;
  mesh: Mesh;
  groupCount: number;
  /** a color per group (atom) */
  colors: DataArray<Uint8Array>;
  /** grid spacing, Å */
  resolution: number;
};

/** A color gradient over the atoms, as the `element-index` theme */
function groupColors(groupCount: number) {
  const image = createDataArray(groupCount, 4, Uint8Array);
  for (let i = 0; i < groupCount; ++i) {
    const t = i / Math.max(1, groupCount - 1);
    image.array[i * 4] = Math.round(255 * t);
    image.array[i * 4 + 1] = Math.round(255 * (1 - Math.abs(2 * t - 1)));
    image.array[i * 4 + 2] = Math.round(255 * (1 - t));
    image.array[i * 4 + 3] = 255;
  }
  return image;
}

/** The structure Gaussian surface mesh, as the CPU mesh visual builds it */
async function gaussianSurfaceMesh(atoms: number) {
  const structure = await globuleStructure(atoms);
  const props = { ...PD.getDefaultValues(GaussianDensityParams), resolution: SurfaceResolution };
  const sizeTheme = getPhysicalSizeTheme(structure);
  const { transform, field, idField, radiusFactor, maxRadius } = await computeStructureGaussianDensity(
    structure,
    sizeTheme,
    props,
  ).run();
  const isoLevel = Math.exp(-props.smoothness) / radiusFactor;
  const mesh = await computeMarchingCubesMesh({ isoLevel, scalarField: field, idField }).run();
  Mesh.transform(mesh, transform);
  const extraRadius = props.radiusOffset * (1 + Math.exp(-props.smoothness));
  mesh.setBoundingSphere(Sphere3D.expand(Sphere3D(), structure.boundary.sphere, maxRadius + extraRadius));
  return { structure, mesh, groupCount: structure.elementCount };
}

const Identity = new Float32Array(Mat4.id);
const Instance = new Float32Array([0]);

/** Reused per graphics context and input, as the representation reuses its color grid texture */
const outputs = new WeakMap<GraphicsContext, Texture>();
function getOutput(graphics: GraphicsContext) {
  let texture = outputs.get(graphics);
  if (!texture) {
    texture = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    outputs.set(graphics, texture);
  }
  return texture;
}

/** The CPU mesh visual: accumulation on the CPU, the grid uploaded into a texture */
async function cpuPath(graphics: GraphicsContext, input: SmoothingInput): Promise<PathRun> {
  const { mesh, groupCount, colors, resolution } = input;
  const timer = new StageTimer(graphics);
  const result = calcMeshColorSmoothing(
    {
      vertexCount: mesh.vertexCount,
      instanceCount: 1,
      groupCount,
      transformBuffer: Identity,
      instanceBuffer: Instance,
      positionBuffer: mesh.vertexBuffer.ref.value,
      groupBuffer: mesh.groupBuffer.ref.value,
      colorData: colors,
      colorType: 'group',
      boundingSphere: mesh.boundingSphere,
      invariantBoundingSphere: mesh.boundingSphere,
      itemSize: 4,
    },
    { resolution, stride: SampleStride },
    graphics,
    getOutput(graphics),
  );
  await timer.stage('smoothing');
  if (result.kind !== 'volume') throw new Error('expected a color grid');
  return { stages: timer.stages, triangles: mesh.triangleCount, bytes: result.texture.getByteCount() };
}

/** The mesh as texture mesh textures (unindexed), as the GPU isosurface outputs them */
type TextureMeshData = { vertexCount: number; position: Texture; group: Texture; colors: Texture };
const textureMeshData = new WeakMap<SmoothingInput, WeakMap<GraphicsContext, TextureMeshData>>();

function getTextureMeshData(graphics: GraphicsContext, input: SmoothingInput): TextureMeshData {
  let perGraphics = textureMeshData.get(input);
  if (!perGraphics) {
    perGraphics = new WeakMap();
    textureMeshData.set(input, perGraphics);
  }
  let data = perGraphics.get(graphics);
  if (!data) {
    const { mesh, colors } = input;
    const positions = mesh.vertexBuffer.ref.value;
    const groups = mesh.groupBuffer.ref.value;
    const indices = mesh.indexBuffer.ref.value;
    const vertexCount = mesh.triangleCount * 3;
    const width = Math.ceil(Math.sqrt(vertexCount));
    const height = Math.ceil(vertexCount / width);
    const position = new Float32Array(width * height * 4);
    const group = new Uint8Array(width * height * 4);
    const rgb = [0, 0, 0];
    for (let i = 0; i < vertexCount; ++i) {
      const v = indices[i];
      for (let c = 0; c < 3; ++c) position[i * 4 + c] = positions[v * 3 + c];
      packIntToRGBArray(groups[v], rgb, 0);
      for (let c = 0; c < 3; ++c) group[i * 4 + c] = rgb[c];
    }
    const { resources } = graphics;
    data = {
      vertexCount,
      position: resources.texture('image-float32', 'rgba', 'float', 'nearest'),
      group: resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest'),
      colors: resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest'),
    };
    data.position.load({ array: position, width, height });
    data.group.load({ array: group, width, height });
    data.colors.load(getDataTextureImage(colors, 4));
    perGraphics.set(graphics, data);
  }
  return data;
}

/** The WebGL texture mesh visual: accumulation with point sprites */
async function gpuPath(graphics: GraphicsContext, input: SmoothingInput): Promise<PathRun> {
  const { mesh, groupCount, resolution } = input;
  const data = getTextureMeshData(graphics, input);
  const timer = new StageTimer(graphics);
  const result = calcTextureMeshColorSmoothing(
    {
      vertexCount: data.vertexCount,
      instanceCount: 1,
      groupCount,
      transformBuffer: Identity,
      instanceBuffer: Instance,
      positionTexture: data.position,
      groupTexture: data.group,
      colorData: data.colors,
      colorType: 'group',
      boundingSphere: mesh.boundingSphere,
      invariantBoundingSphere: mesh.boundingSphere,
    },
    // texture meshes are never indexed, so the stride is tripled
    { resolution, stride: SampleStride * 3 },
    graphics,
    getOutput(graphics),
  );
  await timer.stage('smoothing');
  return { stages: timer.stages, triangles: mesh.triangleCount, bytes: result.texture.getByteCount() };
}

async function webglPath(graphics: GraphicsContext, input: SmoothingInput): Promise<PathRun | string> {
  if (!isWebGLContext(graphics)) return 'needs a WebGL context';
  if (!graphics.capabilities.multipleRenderTargets) return 'needs multiple render targets';
  return gpuPath(graphics, input);
}

const NullTexture = createNullTexture();

/** The surface in the vertex arena, as the WebGPU visual extracts it (GPU density, compute marching cubes) */
type ChunkedData = { surface: ChunkedSurface; vertexCount: number; colors: Texture };
const chunkedData = new WeakMap<SmoothingInput, WeakMap<GraphicsContext, ChunkedData>>();

async function getChunkedData(graphics: GraphicsContext, input: SmoothingInput): Promise<ChunkedData> {
  let perGraphics = chunkedData.get(input);
  if (!perGraphics) {
    perGraphics = new WeakMap();
    chunkedData.set(input, perGraphics);
  }
  let data = perGraphics.get(graphics);
  if (!data) {
    const props = { ...PD.getDefaultValues(GaussianDensityParams), resolution: SurfaceResolution };
    const sizeTheme = getPhysicalSizeTheme(input.structure);
    const density = computeStructureGaussianDensityTexture(input.structure, sizeTheme, props, graphics);
    const isoLevel = Math.exp(-props.smoothness) / density.radiusFactor;
    const result = await computeChunkedIsosurface(
      graphics.compute!,
      { kind: 'texture', texture: density.texture, channel: 'alpha', packedGroup: true },
      density.gridDim,
      isoLevel,
      density.transform,
      false,
    );
    const { drawn } = await result.count;
    const colors = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
    colors.load(getDataTextureImage(input.colors, 4));
    data = { surface: result.surface, vertexCount: drawn, colors };
    perGraphics.set(graphics, data);
  }
  return data;
}

/** The WebGPU texture mesh visual: accumulation with compute kernels, the surface in the vertex arena */
async function webgpuPath(graphics: GraphicsContext, input: SmoothingInput): Promise<PathRun | string> {
  if (!isComputeColorSmoothingSupported(graphics) || !graphics.compute) return 'needs compute shaders';
  const { mesh, groupCount, resolution } = input;
  const data = await getChunkedData(graphics, input);
  const timer = new StageTimer(graphics);
  const result = calcTextureMeshColorSmoothing(
    {
      vertexCount: data.vertexCount,
      instanceCount: 1,
      groupCount,
      transformBuffer: Identity,
      instanceBuffer: Instance,
      // as the values of a texture mesh in the vertex arena
      positionTexture: NullTexture,
      groupTexture: NullTexture,
      chunks: { table: data.surface.table, vertices: data.surface.arena.vertices.ref.value },
      colorData: data.colors,
      colorType: 'group',
      boundingSphere: mesh.boundingSphere,
      invariantBoundingSphere: mesh.boundingSphere,
    },
    // texture meshes are never indexed, so the stride is tripled
    { resolution, stride: SampleStride * 3 },
    graphics,
    getOutput(graphics),
  );
  await timer.stage('smoothing');
  return { stages: timer.stages, triangles: data.vertexCount / 3, bytes: result.texture.getByteCount() };
}

function smoothingCase(atoms: number, resolutionFactor: number): BenchmarkCase<SmoothingInput> {
  const label = atoms >= 1000 ? `${atoms / 1000}k` : `${atoms}`;
  const resolution = SurfaceResolution * resolutionFactor;
  return {
    id: `smoothing/gaussian/globule-${label}-res${resolution}`,
    description:
      `Color smoothing on the Gaussian surface (resolution ${SurfaceResolution} Å) of a protein-like ball of ` +
      `${atoms} atoms, a color per atom, on a ${resolution} Å grid, every ${SampleStride}rd vertex sampled.`,
    setup: async () => {
      const { structure, mesh, groupCount } = await gaussianSurfaceMesh(atoms);
      return { structure, mesh, groupCount, colors: groupColors(groupCount), resolution };
    },
    paths: { cpu: cpuPath, webgl: webglPath, webgpu: webgpuPath },
  };
}

export const ColorSmoothingCases: BenchmarkCase<any>[] = [1000, 10000, 100000].flatMap((atoms) => [
  smoothingCase(atoms, 2),
  smoothingCase(atoms, 1),
]);
