import type { GaussianDensityTextureData } from '../gaussian-density.js';
/**
 * Copyright (c) 2017-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author Michael Krone <michael.krone@uni-tuebingen.de>
 */

import type { PositionData } from '@molstar/core/math/geometry/common';
import { Box3D } from '@molstar/core/math/geometry';
import type { GaussianDensityProps } from '@molstar/core/math/geometry/gaussian-density';
import { OrderedSet } from '@molstar/core/data/int';
import { Vec3, Mat4, Vec2 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { createComputeRenderable, type ComputeRenderable } from '@molstar/graphics/gl/renderable';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { BlendState, type RenderState } from '@molstar/graphics/gl/render-state';
import type { Texture, TextureFilter, TextureFormat, TextureKind, TextureType } from '@molstar/graphics/gl/texture';
import {
  computeGaussianDensity,
  isComputeGaussianDensitySupported,
} from '@molstar/graphics/gl/compute/gaussian-density';
import { ShaderCode } from '@molstar/graphics/gl/shader-code';
import {
  ValueSpec,
  AttributeSpec,
  UniformSpec,
  TextureSpec,
  DefineSpec,
  type Values,
} from '@molstar/graphics/gl/renderable/schema';
import type { Framebuffer } from '@molstar/graphics/gl/framebuffer';
import { isTimingMode } from '@molstar/core/util/debug';

export const GaussianDensitySchema = {
  drawCount: ValueSpec('number'),
  instanceCount: ValueSpec('number'),

  /** quad corners, one instance per atom */
  aMapping: AttributeSpec('float32', 2, 0),
  aRadius: AttributeSpec('float32', 1, 1),
  aPosition: AttributeSpec('float32', 3, 1),
  aGroup: AttributeSpec('float32', 1, 1),

  uCurrentSlice: UniformSpec('f'),
  uCurrentX: UniformSpec('f'),
  uCurrentY: UniformSpec('f'),
  uBboxMin: UniformSpec('v3', 'material'),
  uBboxSize: UniformSpec('v3', 'material'),
  uGridDim: UniformSpec('v3', 'material'),
  uGridTexDim: UniformSpec('v3', 'material'),
  uGridTexScale: UniformSpec('v2', 'material'),
  uAlpha: UniformSpec('f', 'material'),
  uResolution: UniformSpec('f', 'material'),
  uRadiusFactorInv: UniformSpec('f', 'material'),
  tMinDistanceTex: TextureSpec('texture', 'rgba', 'float', 'nearest', 'material', {
    dimension: { define: 'dGridTexType' },
  }),

  dGridTexType: DefineSpec('string', ['2d', '3d']),
  dCalcType: DefineSpec('string', ['density', 'minDistance', 'groupId']),
};
type GaussianDensityValues = Values<typeof GaussianDensitySchema>;
type GaussianDensityRenderable = ComputeRenderable<GaussianDensityValues>;
const GaussianDensityName = 'gaussian-density';
export const GaussianDensityShaderCode = ShaderCode(GaussianDensityName);

/** Whether the GPU gaussian density is available: the splatting program, or compute kernels for the 3d texture */
export function isGaussianDensityTextureSupported(graphics: GraphicsContext) {
  return graphics.isProgramSupported(GaussianDensityShaderCode) || isComputeGaussianDensitySupported(graphics);
}

function getFramebuffer(graphics: GraphicsContext): Framebuffer {
  if (!graphics.namedFramebuffers[GaussianDensityName]) {
    graphics.namedFramebuffers[GaussianDensityName] = graphics.resources.framebuffer();
  }
  return graphics.namedFramebuffers[GaussianDensityName];
}

function getTexture(
  name: string,
  graphics: GraphicsContext,
  kind: TextureKind,
  format: TextureFormat,
  type: TextureType,
  filter: TextureFilter,
): Texture {
  const _name = `${GaussianDensityName}-${name}`;
  if (!graphics.namedTextures[_name]) {
    graphics.namedTextures[_name] = graphics.resources.texture(kind, format, type, filter);
  }
  return graphics.namedTextures[_name];
}

/**
 * `growTexture`: with compute kernels, `oldTexture` is only grown (a texture shared by grids of varying size, read by
 * grid coordinates; `gridTexDim` is then its size)
 */
export function GaussianDensityTexture(
  graphics: GraphicsContext,
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  oldTexture?: Texture,
  growTexture?: boolean,
): GaussianDensityTextureData {
  return graphics.capabilities.texture3d
    ? GaussianDensityTexture3d(graphics, position, box, radius, props, oldTexture, growTexture)
    : GaussianDensityTexture2d(graphics, position, box, radius, false, props, oldTexture);
}

export function GaussianDensityTexture2d(
  graphics: GraphicsContext,
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  powerOfTwo: boolean,
  props: GaussianDensityProps,
  oldTexture?: Texture,
): GaussianDensityTextureData {
  if (isTimingMode) graphics.timer.group('GaussianDensityTexture2d');
  const data = calcGaussianDensityTexture2d(graphics, position, box, radius, powerOfTwo, props, oldTexture);
  if (isTimingMode) graphics.timer.groupEnd('GaussianDensityTexture2d');
  return finalizeGaussianDensityTexture(data);
}

export function GaussianDensityTexture3d(
  graphics: GraphicsContext,
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  oldTexture?: Texture,
  growTexture?: boolean,
): GaussianDensityTextureData {
  if (isTimingMode) graphics.timer.group('GaussianDensityTexture3d');
  const data = calcGaussianDensityTexture3d(graphics, position, box, radius, props, oldTexture, growTexture);
  if (isTimingMode) graphics.timer.groupEnd('GaussianDensityTexture3d');
  return finalizeGaussianDensityTexture(data);
}

function finalizeGaussianDensityTexture({
  texture,
  scale,
  bbox,
  gridDim,
  gridTexDim,
  gridDataDim,
  gridTexScale,
  radiusFactor,
  resolution,
  maxRadius,
}: _GaussianDensityTextureData): GaussianDensityTextureData {
  return {
    transform: getTransform(scale, bbox),
    texture,
    bbox,
    gridDim,
    gridTexDim,
    gridDataDim,
    gridTexScale,
    radiusFactor,
    resolution,
    maxRadius,
  };
}

function getTransform(scale: Vec3, bbox: Box3D) {
  const transform = Mat4.identity();
  Mat4.fromScaling(transform, scale);
  Mat4.setTranslation(transform, bbox.min);
  return transform;
}

//

type _GaussianDensityTextureData = {
  texture: Texture;
  scale: Vec3;
  bbox: Box3D;
  gridDim: Vec3;
  gridTexDim: Vec3;
  gridDataDim: Vec3;
  gridTexScale: Vec2;
  radiusFactor: number;
  resolution: number;
  maxRadius: number;
};

function calcGaussianDensityTexture2d(
  graphics: GraphicsContext,
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  powerOfTwo: boolean,
  props: GaussianDensityProps,
  texture?: Texture,
): _GaussianDensityTextureData {
  // console.log('2d');
  const {
    resources,
    capabilities: { floatBlend, halfFloatRenderTarget },
  } = graphics;
  const { smoothness, resolution } = props;

  const { drawCount, positions, radii, groups, scale, expandedBox, dim, maxRadius } = prepareGaussianDensityData(
    position,
    box,
    radius,
    props,
  );
  const [dx, dy, dz] = dim;
  const { texDimX, texDimY, texCols, powerOfTwoSize } = getTexture2dSize(dim);
  // console.log({ texDimX, texDimY, texCols, powerOfTwoSize, dim });
  const gridTexDim = Vec3.create(texDimX, texDimY, 0);
  const gridTexScale = Vec2.create(texDimX / powerOfTwoSize, texDimY / powerOfTwoSize);
  const radiusFactor = maxRadius * 2;

  const width = powerOfTwo ? powerOfTwoSize : texDimX;
  const height = powerOfTwo ? powerOfTwoSize : texDimY;

  const minDistTex = getTexture('min-dist-2d', graphics, 'image-uint8', 'rgba', 'ubyte', 'nearest');
  minDistTex.define(width, height);

  const renderable = getGaussianDensityRenderable(
    graphics,
    drawCount,
    positions,
    radii,
    groups,
    minDistTex,
    expandedBox,
    dim,
    gridTexDim,
    gridTexScale,
    smoothness,
    resolution,
    radiusFactor,
  );

  //

  const { uCurrentSlice, uCurrentX, uCurrentY } = renderable.values;

  const framebuffer = getFramebuffer(graphics);

  if (!texture)
    texture = halfFloatRenderTarget
      ? resources.texture('image-float16', 'rgba', 'fp16', 'linear')
      : floatBlend
        ? resources.texture('image-float32', 'rgba', 'float', 'linear')
        : resources.texture('image-uint8', 'rgba', 'ubyte', 'linear');
  texture.define(width, height);

  // console.log(renderable)

  function render(fbTex: Texture, clear: boolean) {
    fbTex.attachFramebuffer(framebuffer, 0);
    graphics.beginRenderPass(framebuffer, { label: 'gaussian density' });
    if (clear) {
      graphics.setViewport(0, 0, width, height);
      graphics.setScissor(0, 0, width, height);
      graphics.clear({ color: [0, 0, 0, 0] });
    }
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
      renderable.render();
      ++currCol;
      currX += dx;
    }
    graphics.flush();
  }

  setupDensityRendering(graphics, renderable);
  render(texture, true);

  setupMinDistanceRendering(graphics, renderable);
  render(minDistTex, true);

  setupGroupIdRendering(graphics, renderable);
  render(texture, false);

  // printTextureImage(readTexture(graphics, minDistTex), { scale: 0.75 });

  return {
    texture,
    scale,
    bbox: expandedBox,
    gridDim: dim,
    gridTexDim,
    gridDataDim: dim,
    gridTexScale,
    radiusFactor,
    resolution,
    maxRadius,
  };
}

function calcGaussianDensityTexture3d(
  graphics: GraphicsContext,
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
  texture?: Texture,
  growTexture?: boolean,
): _GaussianDensityTextureData {
  // console.log('3d');
  const {
    resources,
    capabilities: { floatBlend, halfFloatRenderTarget },
  } = graphics;
  const { smoothness, resolution } = props;

  const { drawCount, positions, radii, groups, scale, expandedBox, dim, maxRadius } = prepareGaussianDensityData(
    position,
    box,
    radius,
    props,
  );
  const [dx, dy, dz] = dim;
  const gridTexScale = Vec2.create(1, 1);
  const radiusFactor = maxRadius * 2;

  if (isComputeGaussianDensitySupported(graphics)) {
    if (!texture) texture = resources.texture('volume-float16', 'rgba', 'fp16', 'linear');
    computeGaussianDensity(
      graphics.compute,
      { count: drawCount, positions, radii, groups, maxRadius },
      { box: expandedBox, dim, resolution, smoothness, radiusFactor },
      texture,
      { growTexture },
    );
    return {
      texture,
      scale,
      bbox: expandedBox,
      gridDim: dim,
      gridTexDim: growTexture ? Vec3.create(texture.getWidth(), texture.getHeight(), texture.getDepth()) : dim,
      gridDataDim: dim,
      gridTexScale,
      radiusFactor,
      resolution,
      maxRadius,
    };
  }

  const minDistTex = getTexture('min-dist-3d', graphics, 'volume-uint8', 'rgba', 'ubyte', 'nearest');
  minDistTex.define(dx, dy, dz);

  const renderable = getGaussianDensityRenderable(
    graphics,
    drawCount,
    positions,
    radii,
    groups,
    minDistTex,
    expandedBox,
    dim,
    dim,
    gridTexScale,
    smoothness,
    resolution,
    radiusFactor,
  );

  //

  const { uCurrentSlice } = renderable.values;

  const framebuffer = getFramebuffer(graphics);
  graphics.setViewport(0, 0, dx, dy);
  graphics.setScissor(0, 0, dx, dy);

  if (!texture)
    texture = halfFloatRenderTarget
      ? resources.texture('volume-float16', 'rgba', 'fp16', 'linear')
      : floatBlend
        ? resources.texture('volume-float32', 'rgba', 'float', 'linear')
        : resources.texture('volume-uint8', 'rgba', 'ubyte', 'linear');
  texture.define(dx, dy, dz);

  function render(fbTex: Texture, clear: boolean) {
    for (let i = 0; i < dz; ++i) {
      ValueCell.update(uCurrentSlice, i);
      fbTex.attachFramebuffer(framebuffer, 0, i);
      graphics.beginRenderPass(framebuffer, {
        label: 'gaussian density',
        ...(clear ? { clear: { color: [0, 0, 0, 0] } } : {}),
      });
      renderable.render();
    }
    graphics.flush();
  }

  setupDensityRendering(graphics, renderable);
  render(texture, true);

  setupMinDistanceRendering(graphics, renderable);
  render(minDistTex, true);

  setupGroupIdRendering(graphics, renderable);
  render(texture, false);

  return {
    texture,
    scale,
    bbox: expandedBox,
    gridDim: dim,
    gridTexDim: dim,
    gridDataDim: dim,
    gridTexScale,
    radiusFactor,
    resolution,
    maxRadius,
  };
}

//

function prepareGaussianDensityData(
  position: PositionData,
  box: Box3D,
  radius: (index: number) => number,
  props: GaussianDensityProps,
) {
  const { resolution, radiusOffset } = props;
  const scaleFactor = 1 / resolution;

  const { indices, x, y, z, id } = position;
  const n = OrderedSet.size(indices);

  const positions = new Float32Array(n * 3);
  const radii = new Float32Array(n);
  const groups = new Float32Array(n);

  let maxRadius = 0;

  // a range (the common case) indexed directly, not by `OrderedSet.getAt` per element
  const offset = OrderedSet.isInterval(indices) ? OrderedSet.start(indices) : -1;
  for (let i = 0; i < n; ++i) {
    const j = offset >= 0 ? offset + i : OrderedSet.getAt(indices, i);

    positions[i * 3] = x[j];
    positions[i * 3 + 1] = y[j];
    positions[i * 3 + 2] = z[j];
    const r = radius(j) + radiusOffset;
    if (maxRadius < r) maxRadius = r;
    radii[i] = r;
    groups[i] = id ? id[i] : i;
  }

  const pad = maxRadius * 2 + resolution * 4;
  const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));
  const scaledBox = Box3D.scale(Box3D(), expandedBox, scaleFactor);
  const dim = Box3D.size(Vec3(), scaledBox);
  Vec3.ceil(dim, dim);

  const scale = Vec3.create(resolution, resolution, resolution);

  return { drawCount: n, positions, radii, groups, scale, expandedBox, dim, maxRadius };
}

function getGaussianDensityRenderable(
  graphics: GraphicsContext,
  drawCount: number,
  positions: Float32Array,
  radii: Float32Array,
  groups: Float32Array,
  minDistanceTexture: Texture,
  box: Box3D,
  gridDim: Vec3,
  gridTexDim: Vec3,
  gridTexScale: Vec2,
  smoothness: number,
  resolution: number,
  radiusFactor: number,
) {
  // console.log('radiusFactor', radiusFactor);
  if (graphics.namedComputeRenderables[GaussianDensityName]) {
    const extent = Vec3.sub(Vec3(), box.max, box.min);
    const v = graphics.namedComputeRenderables[GaussianDensityName].values;

    ValueCell.updateIfChanged(v.drawCount, 6);
    ValueCell.updateIfChanged(v.instanceCount, drawCount);

    ValueCell.update(v.aRadius, radii);
    ValueCell.update(v.aPosition, positions);
    ValueCell.update(v.aGroup, groups);

    ValueCell.updateIfChanged(v.uCurrentSlice, 0);
    ValueCell.updateIfChanged(v.uCurrentX, 0);
    ValueCell.updateIfChanged(v.uCurrentY, 0);
    ValueCell.update(v.uBboxMin, box.min);
    ValueCell.update(v.uBboxSize, extent);
    ValueCell.update(v.uGridDim, gridDim);
    ValueCell.update(v.uGridTexDim, gridTexDim);
    ValueCell.update(v.uGridTexScale, gridTexScale);
    ValueCell.updateIfChanged(v.uAlpha, smoothness);
    ValueCell.updateIfChanged(v.uResolution, resolution);
    ValueCell.updateIfChanged(v.uRadiusFactorInv, 1 / radiusFactor);
    ValueCell.update(v.tMinDistanceTex, minDistanceTexture);

    ValueCell.updateIfChanged(v.dGridTexType, minDistanceTexture.getDepth() > 0 ? '3d' : '2d');
    ValueCell.updateIfChanged(v.dCalcType, 'density');

    graphics.namedComputeRenderables[GaussianDensityName].update();
  } else {
    graphics.namedComputeRenderables[GaussianDensityName] = createGaussianDensityRenderable(
      graphics,
      drawCount,
      positions,
      radii,
      groups,
      minDistanceTexture,
      box,
      gridDim,
      gridTexDim,
      gridTexScale,
      smoothness,
      resolution,
      radiusFactor,
    );
  }
  return graphics.namedComputeRenderables[GaussianDensityName];
}

function createGaussianDensityRenderable(
  graphics: GraphicsContext,
  drawCount: number,
  positions: Float32Array,
  radii: Float32Array,
  groups: Float32Array,
  minDistanceTexture: Texture,
  box: Box3D,
  gridDim: Vec3,
  gridTexDim: Vec3,
  gridTexScale: Vec2,
  smoothness: number,
  resolution: number,
  radiusFactor: number,
) {
  const extent = Vec3.sub(Vec3(), box.max, box.min);

  const values: GaussianDensityValues = {
    drawCount: ValueCell.create(6),
    instanceCount: ValueCell.create(drawCount),

    aMapping: ValueCell.create(QuadMapping),
    aRadius: ValueCell.create(radii),
    aPosition: ValueCell.create(positions),
    aGroup: ValueCell.create(groups),

    uCurrentSlice: ValueCell.create(0),
    uCurrentX: ValueCell.create(0),
    uCurrentY: ValueCell.create(0),
    uBboxMin: ValueCell.create(box.min),
    uBboxSize: ValueCell.create(extent),
    uGridDim: ValueCell.create(gridDim),
    uGridTexDim: ValueCell.create(gridTexDim),
    uGridTexScale: ValueCell.create(gridTexScale),
    uAlpha: ValueCell.create(smoothness),
    uResolution: ValueCell.create(resolution),
    uRadiusFactorInv: ValueCell.create(1 / radiusFactor),
    tMinDistanceTex: ValueCell.create(minDistanceTexture),

    dGridTexType: ValueCell.create(minDistanceTexture.getDepth() > 0 ? '3d' : '2d'),
    dCalcType: ValueCell.create('density'),
  };

  const schema = { ...GaussianDensitySchema };
  const renderItem = graphics.createComputeRenderItem('triangles', GaussianDensityShaderCode, schema, values);

  return createComputeRenderable(renderItem, values);
}

/** corners of the two triangles of a quad */
const QuadMapping = new Float32Array([-1, 1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1]);

/** the shader writes 1 - dist so we set blending to MAX */
const MinDistanceState: RenderState = { blend: BlendState.Max, colorWrite: [false, false, false, true] };
const DensityState: RenderState = { blend: BlendState.Additive, colorWrite: [false, false, false, true] };
/** overwrite color, don't change alpha */
const GroupIdState: RenderState = { colorWrite: [true, true, true, false] };

function setupMinDistanceRendering(graphics: GraphicsContext, renderable: GaussianDensityRenderable) {
  ValueCell.update(renderable.values.dCalcType, 'minDistance');
  renderable.update();
  graphics.setRenderState(MinDistanceState);
}

function setupDensityRendering(graphics: GraphicsContext, renderable: GaussianDensityRenderable) {
  ValueCell.update(renderable.values.dCalcType, 'density');
  renderable.update();
  graphics.setRenderState(DensityState);
}

function setupGroupIdRendering(graphics: GraphicsContext, renderable: GaussianDensityRenderable) {
  ValueCell.update(renderable.values.dCalcType, 'groupId');
  renderable.update();
  graphics.setRenderState(GroupIdState);
}

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
