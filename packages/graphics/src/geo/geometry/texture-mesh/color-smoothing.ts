/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { setColorType } from '../color-data.js';
import { setEmissiveType } from '../emissive-data.js';
import { setOverpaintType } from '../overpaint-data.js';
import { setSubstanceType } from '../substance-data.js';
import { setTransparencyType } from '../transparency-data.js';
import { ValueCell } from '@molstar/core/util';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { RenderState } from '@molstar/graphics/gl/render-state';
import { isNullTexture, type Texture } from '@molstar/graphics/gl/texture';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import {
  ValueSpec,
  AttributeSpec,
  UniformSpec,
  TextureSpec,
  type Values,
  DefineSpec,
} from '@molstar/graphics/gl/renderable/schema';
import { QuadSchema, QuadValues } from '@molstar/graphics/gl/compute/util';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Box3D, type Sphere3D } from '@molstar/core/math/geometry';
import type { TextureMeshValues } from '@molstar/graphics/gl/renderable/texture-mesh';
import { TextureMesh } from './texture-mesh.js';
import { isTimingMode } from '@molstar/core/util/debug';
import { computeColorSmoothing, isComputeColorSmoothingSupported } from '@molstar/graphics/gl/compute/color-smoothing';
import type { StorageBuffer } from '@molstar/graphics/gl/compute';
import { type DataArray, getDataTextureImage } from '@molstar/graphics/gl/renderable/util';
import { getRowReachedTest, getSplatInstances } from '../mesh/color-smoothing.js';

export const ColorAccumulateSchema = {
  drawCount: ValueSpec('number'),
  instanceCount: ValueSpec('number'),
  stride: ValueSpec('number'),

  uGroupCount: UniformSpec('i', 'material'),

  aTransform: AttributeSpec('float32', 16, 1),
  aInstance: AttributeSpec('float32', 1, 1),
  aSample: AttributeSpec('float32', 1, 0),

  uGeoTexDim: UniformSpec('v2', 'material'),
  tPosition: TextureSpec('texture', 'rgba', 'float', 'nearest', 'material'),
  tGroup: TextureSpec('texture', 'rgba', 'float', 'nearest', 'material'),

  uColorTexDim: UniformSpec('v2', 'material'),
  tColor: TextureSpec('texture', 'rgba', 'ubyte', 'nearest', 'material'),
  dColorType: DefineSpec('string', ['group', 'groupInstance', 'vertex', 'vertexInstance']),

  uCurrentSlice: UniformSpec('f'),
  uCurrentX: UniformSpec('f'),
  uCurrentY: UniformSpec('f'),
  uBboxMin: UniformSpec('v3', 'material'),
  uBboxSize: UniformSpec('v3', 'material'),
  uResolution: UniformSpec('f', 'material'),
};
export const ColorAccumulateShaderCode = ShaderCode('color-smoothing-accumulate', {
  drawBuffers: 'required',
});
type ColorAccumulateValues = Values<typeof ColorAccumulateSchema>;
const ColorAccumulateName = 'color-accumulate';
const ColorCountName = 'color-count';

interface AccumulateInput {
  /** an upper bound for vertices in the vertex arena (`chunks`) */
  vertexCount: number;
  instanceCount: number;
  groupCount: number;
  transformBuffer: Float32Array;
  instanceBuffer: Float32Array;
  positionTexture: Texture;
  groupTexture: Texture;
  /** the page table and positions of a surface in the vertex arena (WebGPU), instead of the textures; smoothed with
   * compute kernels */
  chunks?: { table: StorageBuffer; vertices: StorageBuffer };
  colorData: Texture;
  colorType: 'group' | 'groupInstance';
}

/** A box holding every untransformed vertex, see `TextureMesh['meta']['vertexBox']` */
function getVertexBox(values: TextureMeshValues) {
  return (values.meta.ref.value as TextureMesh['meta']).vertexBox;
}

/** The vertices in pages of the vertex arena, see `TextureMesh.isoChunks` */
function getChunks(values: TextureMeshValues): AccumulateInput['chunks'] {
  const table = values.sIsoTable.ref.value;
  const vertices = values.sIsoVertices.ref.value;
  return values.dIsoChunks.ref.value && table && vertices ? { table, vertices } : undefined;
}

function getSampleBuffer(sampleCount: number, stride: number) {
  const sampleBuffer = new Float32Array(sampleCount);
  for (let i = 0; i < sampleCount; ++i) {
    sampleBuffer[i] = i * stride;
  }
  return sampleBuffer;
}

function getAccumulateRenderable(
  ctx: GraphicsContext,
  input: AccumulateInput,
  box: Box3D,
  resolution: number,
  stride: number,
): ComputeRenderable<ColorAccumulateValues> {
  if (ctx.namedComputeRenderables[ColorAccumulateName]) {
    const extent = Vec3.sub(Vec3(), box.max, box.min);
    const v = ctx.namedComputeRenderables[ColorAccumulateName].values as ColorAccumulateValues;

    const sampleCount = Math.round(input.vertexCount / stride);
    if (sampleCount > v.drawCount.ref.value || stride !== v.stride.ref.value) {
      ValueCell.update(v.aSample, getSampleBuffer(sampleCount, stride));
    }

    ValueCell.updateIfChanged(v.drawCount, sampleCount);
    ValueCell.updateIfChanged(v.instanceCount, input.instanceCount);
    ValueCell.updateIfChanged(v.stride, stride);

    ValueCell.updateIfChanged(v.uGroupCount, input.groupCount);

    ValueCell.update(v.aTransform, input.transformBuffer);
    ValueCell.update(v.aInstance, input.instanceBuffer);

    ValueCell.update(
      v.uGeoTexDim,
      Vec2.set(v.uGeoTexDim.ref.value, input.positionTexture.getWidth(), input.positionTexture.getHeight()),
    );
    ValueCell.update(v.tPosition, input.positionTexture);
    ValueCell.update(v.tGroup, input.groupTexture);

    ValueCell.update(
      v.uColorTexDim,
      Vec2.set(v.uColorTexDim.ref.value, input.colorData.getWidth(), input.colorData.getHeight()),
    );
    ValueCell.update(v.tColor, input.colorData);
    ValueCell.updateIfChanged(v.dColorType, input.colorType);

    ValueCell.updateIfChanged(v.uCurrentSlice, 0);
    ValueCell.updateIfChanged(v.uCurrentX, 0);
    ValueCell.updateIfChanged(v.uCurrentY, 0);
    ValueCell.update(v.uBboxMin, box.min);
    ValueCell.update(v.uBboxSize, extent);
    ValueCell.updateIfChanged(v.uResolution, resolution);

    ctx.namedComputeRenderables[ColorAccumulateName].update();
  } else {
    ctx.namedComputeRenderables[ColorAccumulateName] = createAccumulateRenderable(ctx, input, box, resolution, stride);
  }
  return ctx.namedComputeRenderables[ColorAccumulateName];
}

function createAccumulateRenderable(
  ctx: GraphicsContext,
  input: AccumulateInput,
  box: Box3D,
  resolution: number,
  stride: number,
) {
  const extent = Vec3.sub(Vec3(), box.max, box.min);
  const sampleCount = Math.round(input.vertexCount / stride);

  const values: ColorAccumulateValues = {
    drawCount: ValueCell.create(sampleCount),
    instanceCount: ValueCell.create(input.instanceCount),
    stride: ValueCell.create(stride),

    uGroupCount: ValueCell.create(input.groupCount),

    aTransform: ValueCell.create(input.transformBuffer),
    aInstance: ValueCell.create(input.instanceBuffer),
    aSample: ValueCell.create(getSampleBuffer(sampleCount, stride)),

    uGeoTexDim: ValueCell.create(Vec2.create(input.positionTexture.getWidth(), input.positionTexture.getHeight())),
    tPosition: ValueCell.create(input.positionTexture),
    tGroup: ValueCell.create(input.groupTexture),

    uColorTexDim: ValueCell.create(Vec2.create(input.colorData.getWidth(), input.colorData.getHeight())),
    tColor: ValueCell.create(input.colorData),
    dColorType: ValueCell.create(input.colorType),

    uCurrentSlice: ValueCell.create(0),
    uCurrentX: ValueCell.create(0),
    uCurrentY: ValueCell.create(0),
    uBboxMin: ValueCell.create(box.min),
    uBboxSize: ValueCell.create(extent),
    uResolution: ValueCell.create(resolution),
  };

  const schema = { ...ColorAccumulateSchema };
  const renderItem = ctx.createComputeRenderItem('points', ColorAccumulateShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

export const ColorNormalizeSchema = {
  ...QuadSchema,

  tColor: TextureSpec('texture', 'rgba', 'float', 'nearest'),
  tCount: TextureSpec('texture', 'r', 'float', 'nearest'),
  uTexSize: UniformSpec('v2'),
};
export const ColorNormalizeShaderCode = ShaderCode('normalize');
type ColorNormalizeValues = Values<typeof ColorNormalizeSchema>;
const ColorNormalizeName = 'color-normalize';

function getNormalizeRenderable(
  ctx: GraphicsContext,
  color: Texture,
  count: Texture,
): ComputeRenderable<ColorNormalizeValues> {
  if (ctx.namedComputeRenderables[ColorNormalizeName]) {
    const v = ctx.namedComputeRenderables[ColorNormalizeName].values as ColorNormalizeValues;

    ValueCell.update(v.tColor, color);
    ValueCell.update(v.tCount, count);
    ValueCell.update(v.uTexSize, Vec2.set(v.uTexSize.ref.value, color.getWidth(), color.getHeight()));

    ctx.namedComputeRenderables[ColorNormalizeName].update();
  } else {
    ctx.namedComputeRenderables[ColorNormalizeName] = createColorNormalizeRenderable(ctx, color, count);
  }
  return ctx.namedComputeRenderables[ColorNormalizeName];
}

function createColorNormalizeRenderable(ctx: GraphicsContext, color: Texture, count: Texture) {
  const values: ColorNormalizeValues = {
    ...QuadValues,
    tColor: ValueCell.create(color),
    tCount: ValueCell.create(count),
    uTexSize: ValueCell.create(Vec2.create(color.getWidth(), color.getHeight())),
  };

  const schema = { ...ColorNormalizeSchema };
  const renderItem = ctx.createComputeRenderItem('triangles', ColorNormalizeShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

//

function getTexture2dSize(gridDim: Vec3) {
  const area = gridDim[0] * gridDim[1] * gridDim[2];
  const squareDim = Math.sqrt(area);
  const powerOfTwoSize = Math.pow(2, Math.ceil(Math.log(squareDim) / Math.log(2)));

  let texDimX = 0;
  let texDimY = gridDim[1];
  let texRows = 1;
  let texCols = gridDim[2];
  if (powerOfTwoSize < gridDim[0] * gridDim[2]) {
    texCols = Math.floor(powerOfTwoSize / gridDim[0]);
    texRows = Math.ceil(gridDim[2] / texCols);
    texDimX = texCols * gridDim[0];
    texDimY *= texRows;
  } else {
    texDimX = gridDim[0] * gridDim[2];
  }
  // console.log(texDimX, texDimY, texDimY < powerOfTwoSize ? powerOfTwoSize : powerOfTwoSize * 2);
  return {
    texDimX,
    texDimY,
    texRows,
    texCols,
    powerOfTwoSize: texDimY < powerOfTwoSize ? powerOfTwoSize : powerOfTwoSize * 2,
  };
}

interface ColorSmoothingInput extends AccumulateInput {
  boundingSphere: Sphere3D;
  invariantBoundingSphere: Sphere3D;
  /** the layer as loaded into `colorData`, with `itemSize` values per element */
  layer?: { image: DataArray<Uint8Array>; itemSize: 4 | 1 };
  /** a box holding every untransformed vertex (`TextureMesh['meta']['vertexBox']`) */
  vertexBox?: Box3D;
}

/** The instances of `input` to draw into the grid: `transformBuffer`, `instanceBuffer` and `instanceCount` */
type SplatInstances = Pick<AccumulateInput, 'instanceCount' | 'transformBuffer' | 'instanceBuffer'>;

/**
 * The instances to draw into the grid (perf-audit 13/P1, see `getSplatInstances` of the mesh path): those whose layer
 * row is non-zero and those whose cells overlap a reached one's, in their order; every instance without `layer` or
 * `vertexBox`, and for non-instance types. The grid stays the one of every instance.
 *
 * Exact as on the mesh path: the cells of reached instances get every contribution, in the former order (instanced
 * draws are in instance order, the compute path's sums are fixed-point); the other cells get none and are 0, as
 * before (with the zero layer: 0 / weight; now 0 / 0, which the compute normalization writes as 0, and the point
 * sprite path is left out only where its normalization does, see `probeZeroWeightNormalization`). Single-channel layers (transparency, emissive) are exact in red,
 * the channel the shaders read: their `r` texture reads alpha 1, so the grid's alpha, which no shader reads, was 1
 * wherever any instance splatted and is now 0 where only left-out instances did. The vertices read groups from a
 * texture, which the CPU does not check: the groups of gaussian surfaces are in range or -1 (a cell whose voxel has no
 * atom, packed 0), which reads the last item of the row before; so that item counts for a row, and row 0 is always
 * reached.
 */
export function getTextureMeshSplatInstances(
  input: ColorSmoothingInput,
  gridDim: Vec3,
  min: Vec3,
  scaleFactor: number,
): SplatInstances & { splatCount: number } {
  const { instanceCount, transformBuffer, instanceBuffer, layer, vertexBox } = input;
  const all = { instanceCount, transformBuffer, instanceBuffer, splatCount: instanceCount };
  if (!input.colorType.endsWith('Instance') || !layer || !vertexBox) return all;

  const isRowReached = getRowReachedTest(layer.image.array, input.groupCount, layer.itemSize, 1);
  const { splat } = getSplatInstances(
    transformBuffer,
    vertexBox,
    instanceCount,
    true,
    (i) => isRowReached(instanceBuffer[i]),
    gridDim,
    min,
    scaleFactor,
    false,
  );
  let splatCount = 0;
  for (let i = 0; i < instanceCount; ++i) splatCount += splat[i];
  // none (not expected: row 0 is always reached, and the grid holds every instance) would bind empty buffers
  if (splatCount === instanceCount || splatCount === 0) return all;

  const transforms = new Float32Array(splatCount * 16);
  const instances = new Float32Array(splatCount);
  for (let i = 0, k = 0; i < instanceCount; ++i) {
    if (!splat[i]) continue;
    transforms.set(transformBuffer.subarray(i * 16, i * 16 + 16), k * 16);
    instances[k] = instanceBuffer[i];
    k += 1;
  }
  return { instanceCount: splatCount, transformBuffer: transforms, instanceBuffer: instances, splatCount };
}

export type ColorSmoothingOptions = {
  resolution: number;
  stride: number;
};

/** The framebuffer and the accumulation targets (weighted colors, weights) of the point sprite path */
function getAccumulateTargets(graphics: GraphicsContext) {
  const {
    resources,
    capabilities: { halfFloatRenderTarget, mixedAttachmentFormats },
  } = graphics;
  if (!graphics.namedFramebuffers[ColorAccumulateName]) {
    graphics.namedFramebuffers[ColorAccumulateName] = resources.framebuffer();
  }
  const framebuffer = graphics.namedFramebuffers[ColorAccumulateName];

  if (mixedAttachmentFormats) {
    if (!graphics.namedTextures[ColorAccumulateName]) {
      graphics.namedTextures[ColorAccumulateName] = halfFloatRenderTarget
        ? resources.texture('image-float16', 'rgba', 'fp16', 'nearest')
        : resources.texture('image-float32', 'rgba', 'float', 'nearest');
    }

    if (!graphics.namedTextures[ColorCountName]) {
      graphics.namedTextures[ColorCountName] = resources.texture('image-float32', 'r', 'float', 'nearest');
    }
  } else {
    // webgl1 requires consistent bit plane counts
    // this is quite wasteful but good enough for medium size meshes

    if (!graphics.namedTextures[ColorAccumulateName]) {
      graphics.namedTextures[ColorAccumulateName] = resources.texture('image-float32', 'rgba', 'float', 'nearest');
    }

    if (!graphics.namedTextures[ColorCountName]) {
      graphics.namedTextures[ColorCountName] = resources.texture('image-float32', 'rgba', 'float', 'nearest');
    }
  }

  const accumulateTexture = graphics.namedTextures[ColorAccumulateName];
  const countTexture = graphics.namedTextures[ColorCountName];
  return { framebuffer, accumulateTexture, countTexture };
}

/** Per context: whether the point sprite path writes cells without weight as 0, once known */
const ZeroWeightProbes = new WeakMap<GraphicsContext, { result?: boolean; promise: Promise<boolean> }>();

/**
 * Whether the normalization of the point sprite path (WebGL, and WebGPU without compute kernels) writes a cell without
 * weight, 0 / 0, as 0 into its unorm texture: GLSL and the GL specs leave NaN to unorm undefined (D3D, Metal and Vulkan
 * give 0). Leaving out instances (`getTextureMeshSplatInstances`) turns cells their vertices read from 0 / weight into
 * 0 / 0, so it is only exact where this holds. Probed once per context, on a 1×1 target read back.
 */
export function probeZeroWeightNormalization(graphics: GraphicsContext): Promise<boolean> {
  let probe = ZeroWeightProbes.get(graphics);
  if (!probe) {
    const promise = runZeroWeightProbe(graphics).catch(() => false);
    const entry: { result?: boolean; promise: Promise<boolean> } = { promise };
    promise.then((result) => {
      entry.result = result;
    });
    ZeroWeightProbes.set(graphics, entry);
    probe = entry;
  }
  return probe.promise;
}

/** `probeZeroWeightNormalization` once its result is known; false (and the probe started) before */
function isZeroWeightNormalizedToZero(graphics: GraphicsContext) {
  const probe = ZeroWeightProbes.get(graphics);
  if (!probe) void probeZeroWeightNormalization(graphics);
  return probe?.result === true;
}

async function runZeroWeightProbe(graphics: GraphicsContext) {
  const { framebuffer, accumulateTexture, countTexture } = getAccumulateTargets(graphics);
  accumulateTexture.define(1, 1);
  countTexture.define(1, 1);
  accumulateTexture.attachFramebuffer(framebuffer, 0);
  countTexture.attachFramebuffer(framebuffer, 1);
  framebuffer.setDrawBuffers(2);
  graphics.beginRenderPass(framebuffer, { label: 'color smoothing probe clear' });
  graphics.setViewport(0, 0, 1, 1);
  graphics.setScissor(0, 0, 1, 1);
  graphics.clear({ color: [0, 0, 0, 0] });
  accumulateTexture.detachFramebuffer(framebuffer, 0);
  countTexture.detachFramebuffer(framebuffer, 1);
  framebuffer.setDrawBuffers(1);

  // as the normalization below: additive into a cleared rgba8 target
  const target = graphics.resources.texture('image-uint8', 'rgba', 'ubyte', 'nearest');
  try {
    target.define(1, 1);
    const normalizeRenderable = getNormalizeRenderable(graphics, accumulateTexture, countTexture);
    graphics.setRenderState(RenderState.Additive);
    target.attachFramebuffer(framebuffer, 0);
    graphics.beginRenderPass(framebuffer, { label: 'color smoothing probe' });
    graphics.setViewport(0, 0, 1, 1);
    graphics.setScissor(0, 0, 1, 1);
    graphics.clear({ color: [0, 0, 0, 0] });
    normalizeRenderable.render();
    target.detachFramebuffer(framebuffer, 0);
    const { array } = await graphics.readTextureData(target);
    for (let i = 0; i < 4; ++i) if (array[i] !== 0) return false;
    return true;
  } finally {
    target.destroy();
  }
}

/** Grid of the smoothing: the box expanded by `1 + resolution`, two extra cells per axis */
export function getGridDim(box: Box3D, resolution: number) {
  const pad = 1 + resolution;
  const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
  const gridDim = Box3D.size(Vec3(), Box3D.scale(Box3D(), expandedBox, 1 / resolution));
  Vec3.ceil(gridDim, gridDim);
  return Vec3.add(gridDim, gridDim, Vec3.create(2, 2, 2));
}

/**
 * `resolution`, coarsened in steps of 10% until the 2d atlas of the grid (largest side by `getTextureSide`) fits
 * textures of `maxTextureSize` (larger grids would fail to allocate; WebGPU's limit is often half of WebGL's)
 */
export function getFittingResolution(
  box: Box3D,
  resolution: number,
  maxTextureSize: number,
  getTextureSide: (gridDim: Vec3) => number,
) {
  while (getTextureSide(getGridDim(box, resolution)) > maxTextureSize) resolution *= 1.1;
  return resolution;
}

export function calcTextureMeshColorSmoothing(
  input: ColorSmoothingInput,
  options: ColorSmoothingOptions,
  graphics: GraphicsContext,
  texture?: Texture,
) {
  // a texture of another context (e.g. before a backend switch) is not reused
  if (texture && !graphics.resources.ownsTexture(texture)) texture = undefined;
  const { chunks } = input;
  const useCompute = !!chunks && isComputeColorSmoothingSupported(graphics);
  if (!useCompute && !graphics.capabilities.multipleRenderTargets) throw new Error('need multiple render targets');

  if (isTimingMode) graphics.timer.group('calcTextureMeshColorSmoothing');
  const { resources } = graphics;
  const { stride } = options;

  const isInstanceType = input.colorType.endsWith('Instance');
  const box = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);
  const resolution = getFittingResolution(box, options.resolution, resources.limits.maxTextureSize, (gridDim) => {
    const { texDimX, texDimY } = getTexture2dSize(gridDim);
    return Math.max(texDimX, texDimY);
  });
  const pad = 1 + resolution;
  const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));

  if (!isInstanceType) {
    input.instanceCount = 1;
    input.instanceBuffer = new Float32Array([0]);
    input.transformBuffer = new Float32Array(Mat4.id);
  }

  const scaleFactor = 1 / resolution;
  const gridDim = getGridDim(box, resolution);
  const { min } = expandedBox;
  // the grid spans every instance, also where the layer is zero: the shaders read it without a bounds check
  // (`texture3dFrom2dLinear` in WebGL and WebGPU), so a point beyond a smaller grid would read another slice; the point
  // sprite path leaves out instances only where 0 / 0 is written as 0
  const leaveOut = useCompute || isZeroWeightNormalizedToZero(graphics);
  const subset = getTextureMeshSplatInstances(
    leaveOut ? input : { ...input, layer: undefined },
    gridDim,
    min,
    scaleFactor,
  );
  const { instanceCount, transformBuffer, instanceBuffer } = subset;
  const samples = { ...input, instanceCount, transformBuffer, instanceBuffer };

  const [dx, dy, dz] = gridDim;
  const { texDimX: width, texDimY: height, texCols } = getTexture2dSize(gridDim);
  // console.log({ width, height, texCols, gridDim, resolution });
  const gridTransform = Vec4.create(min[0], min[1], min[2], scaleFactor);
  const type = isInstanceType ? 'volumeInstance' : 'volume';

  if (useCompute) {
    if (!texture || isNullTexture(texture)) {
      texture = resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
    }
    const bboxSize = Box3D.size(Vec3(), expandedBox);
    computeColorSmoothing(
      graphics.compute,
      { ...samples, chunks, stride },
      { gridDim, texCols, width, height, bboxMin: min, bboxSize, resolution },
      texture,
    );
    if (isTimingMode) graphics.timer.groupEnd('calcTextureMeshColorSmoothing');
    return { texture, gridDim, gridTexDim: Vec2.create(width, height), gridTransform, type };
  }

  const { framebuffer, accumulateTexture, countTexture } = getAccumulateTargets(graphics);
  accumulateTexture.define(width, height);
  countTexture.define(width, height);

  accumulateTexture.attachFramebuffer(framebuffer, 0);
  countTexture.attachFramebuffer(framebuffer, 1);

  const accumulateRenderable = getAccumulateRenderable(graphics, samples, expandedBox, resolution, stride);

  framebuffer.setDrawBuffers(2);
  graphics.beginRenderPass(framebuffer, { label: 'color accumulate' });

  const { uCurrentSlice, uCurrentX, uCurrentY } = accumulateRenderable.values;

  graphics.setRenderState(RenderState.Additive);
  graphics.setViewport(0, 0, width, height);
  graphics.setScissor(0, 0, width, height);
  graphics.clear({ color: [0, 0, 0, 0] });
  ValueCell.update(uCurrentY, 0);
  let currCol = 0;
  let currY = 0;
  let currX = 0;
  for (let i = 0; i < dz; ++i) {
    if (currCol >= texCols) {
      currCol -= texCols;
      currY += dy;
      currX = 0;
      ValueCell.update(uCurrentY, currY);
    }
    // console.log({ i, currX, currY });
    ValueCell.update(uCurrentX, currX);
    ValueCell.update(uCurrentSlice, i);
    // tiles at texel rows from the bottom (`uCurrentY`), rectangles have a top-left origin
    const tileY = height - currY - dy;
    graphics.setViewport(currX, tileY, dx, dy);
    graphics.setScissor(currX, tileY, dx, dy);
    accumulateRenderable.render();
    ++currCol;
    currX += dx;
  }

  accumulateTexture.detachFramebuffer(framebuffer, 0);
  countTexture.detachFramebuffer(framebuffer, 1);
  framebuffer.setDrawBuffers(1);

  // const accImage = new Float32Array(width * height * 4);
  // accumulateTexture.attachFramebuffer(framebuffer, 0);
  // webgl.readPixels(0, 0, width, height, accImage);
  // console.log(accImage);
  // printTextureImage({ array: accImage, width, height }, { scale: 1 });

  // const cntImage = new Float32Array(width * height * 4);
  // countTexture.attachFramebuffer(framebuffer, 0);
  // webgl.readPixels(0, 0, width, height, cntImage);
  // console.log(cntImage);
  // printTextureImage({ array: cntImage, width, height }, { scale: 1 });

  // normalize

  if (!texture || isNullTexture(texture)) {
    texture = resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
  }
  texture.define(width, height);

  const normalizeRenderable = getNormalizeRenderable(graphics, accumulateTexture, countTexture);

  graphics.setRenderState(RenderState.Additive);
  texture.attachFramebuffer(framebuffer, 0);
  graphics.beginRenderPass(framebuffer, { label: 'color normalize' });
  graphics.setViewport(0, 0, width, height);
  graphics.setScissor(0, 0, width, height);
  graphics.clear({ color: [0, 0, 0, 0] });
  normalizeRenderable.render();

  // const normImage = new Uint8Array(width * height * 4);
  // texture.attachFramebuffer(framebuffer, 0);
  // webgl.readPixels(0, 0, width, height, normImage);
  // console.log(normImage);
  // printTextureImage({ array: normImage, width, height }, { scale: 1 });

  if (isTimingMode) graphics.timer.groupEnd('calcTextureMeshColorSmoothing');

  // printTextureImage(readTexture(webgl, texture), { scale: 0.75, id: `${texture.id}` });

  return { texture, gridDim, gridTexDim: Vec2.create(width, height), gridTransform, type };
}

//

const ColorSmoothingColorName = 'color-smoothing-color';
const ColorSmoothingRgbaName = 'color-smoothing-rgba';
const ColorSmoothingAlphaName = 'color-smoothing-alpha';

function isSupportedColorType(x: string): x is 'group' | 'groupInstance' {
  return x === 'group' || x === 'groupInstance';
}

/**
 * Runs `run` again when the vertices of the texture mesh changed after it ran (extracted again, or more than its upper
 * bound, see `computeTextureMeshIsosurface`): smoothing samples the vertices, which were incomplete before
 */
function setReextractionHandler(values: TextureMeshValues, key: string, run: () => void) {
  TextureMesh.setReextractionHandler(values.meta.ref.value as TextureMesh['meta'], key, run);
}

export function applyTextureMeshColorSmoothing(
  values: TextureMeshValues,
  options: ColorSmoothingOptions,
  graphics: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedColorType(values.colorType.ref.value)) return;
  const type = values.colorType.ref.value;

  const run = () => {
    const tripled = { ...options, stride: options.stride * 3 }; // triple because TextureMesh is never indexed (no elements buffer)

    if (!graphics.namedTextures[ColorSmoothingColorName]) {
      graphics.namedTextures[ColorSmoothingColorName] = graphics.resources.texture(
        'image-uint8',
        'rgba',
        'ubyte',
        'nearest',
      );
    }
    const colorData = graphics.namedTextures[ColorSmoothingColorName];
    colorData.load(getDataTextureImage(values.tColor.ref.value, 4));

    const smoothingData = calcTextureMeshColorSmoothing(
      {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        chunks: getChunks(values),
        colorData,
        colorType: type,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
        // no layer: base colors (alpha 255) reach every instance
      },
      tripled,
      graphics,
      colorTexture,
    );

    // reused when run again
    colorTexture = smoothingData.texture;
    setColorType(values, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.gridTransform);
  };
  run();
  // unless no longer smoothed
  setReextractionHandler(values, 'color', () => {
    if (values.colorType.ref.value.startsWith('volume')) run();
  });
}

function isSupportedOverpaintType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyTextureMeshOverpaintSmoothing(
  values: TextureMeshValues,
  options: ColorSmoothingOptions,
  graphics: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedOverpaintType(values.overpaintType.ref.value)) return;
  const type = values.overpaintType.ref.value;

  const run = () => {
    const tripled = { ...options, stride: options.stride * 3 }; // triple because TextureMesh is never indexed (no elements buffer)

    if (!graphics.namedTextures[ColorSmoothingRgbaName]) {
      graphics.namedTextures[ColorSmoothingRgbaName] = graphics.resources.texture(
        'image-uint8',
        'rgba',
        'ubyte',
        'nearest',
      );
    }
    const colorData = graphics.namedTextures[ColorSmoothingRgbaName];
    colorData.load(getDataTextureImage(values.tOverpaint.ref.value, 4));

    const smoothingData = calcTextureMeshColorSmoothing(
      {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        chunks: getChunks(values),
        colorData,
        colorType: type,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
        layer: { image: values.tOverpaint.ref.value, itemSize: 4 },
        vertexBox: getVertexBox(values),
      },
      tripled,
      graphics,
      colorTexture,
    );

    // reused when run again
    colorTexture = smoothingData.texture;
    setOverpaintType(values, smoothingData.type);
    ValueCell.update(values.tOverpaintGrid, smoothingData.texture);
    ValueCell.update(values.uOverpaintGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uOverpaintGridDim, smoothingData.gridDim);
    ValueCell.update(values.uOverpaintGridTransform, smoothingData.gridTransform);
  };
  run();
  // unless no longer smoothed
  setReextractionHandler(values, 'overpaint', () => {
    if (values.overpaintType.ref.value.startsWith('volume')) run();
  });
}

function isSupportedTransparencyType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyTextureMeshTransparencySmoothing(
  values: TextureMeshValues,
  options: ColorSmoothingOptions,
  graphics: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedTransparencyType(values.transparencyType.ref.value)) return;
  const type = values.transparencyType.ref.value;

  const run = () => {
    const tripled = { ...options, stride: options.stride * 3 }; // triple because TextureMesh is never indexed (no elements buffer)

    if (!graphics.namedTextures[ColorSmoothingAlphaName]) {
      graphics.namedTextures[ColorSmoothingAlphaName] = graphics.resources.texture(
        'image-uint8',
        'r',
        'ubyte',
        'nearest',
      );
    }
    const colorData = graphics.namedTextures[ColorSmoothingAlphaName];
    colorData.load(getDataTextureImage(values.tTransparency.ref.value, 1));

    const smoothingData = calcTextureMeshColorSmoothing(
      {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        chunks: getChunks(values),
        colorData,
        colorType: type,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
        layer: { image: values.tTransparency.ref.value, itemSize: 1 },
        vertexBox: getVertexBox(values),
      },
      tripled,
      graphics,
      colorTexture,
    );

    // reused when run again
    colorTexture = smoothingData.texture;
    setTransparencyType(values, smoothingData.type);
    ValueCell.update(values.tTransparencyGrid, smoothingData.texture);
    ValueCell.update(values.uTransparencyGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uTransparencyGridDim, smoothingData.gridDim);
    ValueCell.update(values.uTransparencyGridTransform, smoothingData.gridTransform);
  };
  run();
  // unless no longer smoothed
  setReextractionHandler(values, 'transparency', () => {
    if (values.transparencyType.ref.value.startsWith('volume')) run();
  });
}

function isSupportedEmissiveType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyTextureMeshEmissiveSmoothing(
  values: TextureMeshValues,
  options: ColorSmoothingOptions,
  graphics: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedEmissiveType(values.emissiveType.ref.value)) return;
  const type = values.emissiveType.ref.value;

  const run = () => {
    const tripled = { ...options, stride: options.stride * 3 }; // triple because TextureMesh is never indexed (no elements buffer)

    if (!graphics.namedTextures[ColorSmoothingAlphaName]) {
      graphics.namedTextures[ColorSmoothingAlphaName] = graphics.resources.texture(
        'image-uint8',
        'r',
        'ubyte',
        'nearest',
      );
    }
    const colorData = graphics.namedTextures[ColorSmoothingAlphaName];
    colorData.load(getDataTextureImage(values.tEmissive.ref.value, 1));

    const smoothingData = calcTextureMeshColorSmoothing(
      {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        chunks: getChunks(values),
        colorData,
        colorType: type,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
        layer: { image: values.tEmissive.ref.value, itemSize: 1 },
        vertexBox: getVertexBox(values),
      },
      tripled,
      graphics,
      colorTexture,
    );

    // reused when run again
    colorTexture = smoothingData.texture;
    setEmissiveType(values, smoothingData.type);
    ValueCell.update(values.tEmissiveGrid, smoothingData.texture);
    ValueCell.update(values.uEmissiveGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uEmissiveGridDim, smoothingData.gridDim);
    ValueCell.update(values.uEmissiveGridTransform, smoothingData.gridTransform);
  };
  run();
  // unless no longer smoothed
  setReextractionHandler(values, 'emissive', () => {
    if (values.emissiveType.ref.value.startsWith('volume')) run();
  });
}

function isSupportedSubstanceType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyTextureMeshSubstanceSmoothing(
  values: TextureMeshValues,
  options: ColorSmoothingOptions,
  graphics: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedSubstanceType(values.substanceType.ref.value)) return;
  const type = values.substanceType.ref.value;

  const run = () => {
    const tripled = { ...options, stride: options.stride * 3 }; // triple because TextureMesh is never indexed (no elements buffer)

    if (!graphics.namedTextures[ColorSmoothingRgbaName]) {
      graphics.namedTextures[ColorSmoothingRgbaName] = graphics.resources.texture(
        'image-uint8',
        'rgba',
        'ubyte',
        'nearest',
      );
    }
    const colorData = graphics.namedTextures[ColorSmoothingRgbaName];
    colorData.load(getDataTextureImage(values.tSubstance.ref.value, 4));

    const smoothingData = calcTextureMeshColorSmoothing(
      {
        vertexCount: values.uVertexCount.ref.value,
        instanceCount: values.uInstanceCount.ref.value,
        groupCount: values.uGroupCount.ref.value,
        transformBuffer: values.aTransform.ref.value,
        instanceBuffer: values.aInstance.ref.value,
        positionTexture: values.tPosition.ref.value,
        groupTexture: values.tGroup.ref.value,
        chunks: getChunks(values),
        colorData,
        colorType: type,
        boundingSphere: values.boundingSphere.ref.value,
        invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
        layer: { image: values.tSubstance.ref.value, itemSize: 4 },
        vertexBox: getVertexBox(values),
      },
      tripled,
      graphics,
      colorTexture,
    );

    // reused when run again
    colorTexture = smoothingData.texture;
    setSubstanceType(values, smoothingData.type);
    ValueCell.update(values.tSubstanceGrid, smoothingData.texture);
    ValueCell.update(values.uSubstanceGridTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uSubstanceGridDim, smoothingData.gridDim);
    ValueCell.update(values.uSubstanceGridTransform, smoothingData.gridTransform);
  };
  run();
  // unless no longer smoothed
  setReextractionHandler(values, 'substance', () => {
    if (values.substanceType.ref.value.startsWith('volume')) run();
  });
}
