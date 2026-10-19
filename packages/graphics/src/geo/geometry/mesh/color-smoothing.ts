/**
 * Copyright (c) 2021-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { setColorType } from '../color-data.js';
import { setEmissiveType } from '../emissive-data.js';
import { setOverpaintType } from '../overpaint-data.js';
import { setSubstanceType } from '../substance-data.js';
import { setTransparencyType } from '../transparency-data.js';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import { createTextureImage, type TextureImage } from '@molstar/graphics/gl/renderable/util';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { Texture } from '@molstar/graphics/gl/texture';
import { Box3D, type Sphere3D } from '@molstar/core/math/geometry';
import { lerp } from '@molstar/core/math/interpolate';
import { Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { getVolumeTexture2dLayout } from '@molstar/graphics/repr/volume/util';
import { getFittingResolution } from '../texture-mesh/color-smoothing.js';
import { ValueCell } from '@molstar/core/util';

interface ColorSmoothingInput {
  vertexCount: number;
  instanceCount: number;
  groupCount: number;
  transformBuffer: Float32Array;
  instanceBuffer: Float32Array;
  positionBuffer: Float32Array;
  groupBuffer: Float32Array;
  colorData: TextureImage<Uint8Array>;
  colorType: 'group' | 'groupInstance';
  boundingSphere: Sphere3D;
  invariantBoundingSphere: Sphere3D;
  itemSize: 4 | 1;
}

export type ColorSmoothingOptions = {
  resolution: number;
  stride: number;
};

export function calcMeshColorSmoothing(
  input: ColorSmoothingInput,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  texture?: Texture,
) {
  // a texture of another context (e.g. before a backend switch) is not reused
  if (texture && graphics && !graphics.resources.ownsTexture(texture)) texture = undefined;
  const { colorType, vertexCount, groupCount, positionBuffer, instanceBuffer, transformBuffer, groupBuffer, itemSize } =
    input;
  const { stride } = options;

  const isInstanceType = colorType.endsWith('Instance');
  const box = Box3D.fromSphere3D(Box3D(), isInstanceType ? input.boundingSphere : input.invariantBoundingSphere);
  // the grid is uploaded into a texture with a graphics context
  const resolution = graphics
    ? getFittingResolution(box, options.resolution, graphics.resources.limits.maxTextureSize, (gridDim) => {
        const { width, height } = getVolumeTexture2dLayout(gridDim);
        return Math.max(width, height);
      })
    : options.resolution;
  const pad = 1 + resolution;
  const expandedBox = Box3D.expand(Box3D(), box, Vec3.create(pad, pad, pad));

  const scaleFactor = 1 / resolution;
  const scaledBox = Box3D.scale(Box3D(), expandedBox, scaleFactor);
  const gridDim = Box3D.size(Vec3(), scaledBox);
  Vec3.ceil(gridDim, gridDim);
  Vec3.add(gridDim, gridDim, Vec3.create(2, 2, 2));
  const { min } = expandedBox;

  const { width, height } = getVolumeTexture2dLayout(gridDim);

  const data = new Float32Array(width * height * itemSize);
  const count = new Float32Array(width * height);

  const grid = new Uint8Array(width * height * itemSize);
  const textureImage: TextureImage<Uint8Array> = { array: grid, width, height, filter: 'linear' };

  const instanceCount = isInstanceType ? input.instanceCount : 1;
  const colors = input.colorData.array;

  const zOff = getSliceOffsets(gridDim, width);

  const p = 2;
  const p2 = p * p;
  const [dimX, dimY, dimZ] = gridDim;
  const v = Vec3();

  for (let i = 0; i < instanceCount; ++i) {
    // - use reordered index for access from GPU
    // - use serial index for access from CPU
    const instanceIndex = graphics && isInstanceType ? instanceBuffer[i] : i;
    for (let j = 0; j < vertexCount; j += stride) {
      Vec3.fromArray(v, positionBuffer, j * 3);
      if (isInstanceType) Vec3.transformMat4Offset(v, v, transformBuffer, 0, 0, i * 16);
      Vec3.sub(v, v, min);
      Vec3.scale(v, v, scaleFactor);
      const vx = v[0];
      const vy = v[1];
      const vz = v[2];

      // vertex mapped to grid
      const x = Math.floor(vx);
      const y = Math.floor(vy);
      const z = Math.floor(vz);

      // group colors
      const ci = (instanceIndex * groupCount + groupBuffer[j]) * itemSize;
      const c0 = colors[ci];
      const c1 = itemSize === 4 ? colors[ci + 1] : 0;
      const c2 = itemSize === 4 ? colors[ci + 2] : 0;
      const c3 = itemSize === 4 ? colors[ci + 3] : 0;

      // Extents of grid to consider for this atom
      const begX = Math.max(0, x - p);
      const begY = Math.max(0, y - p);
      const begZ = Math.max(0, z - p);

      // Add two to these points:
      // - x, y, z are floor'd values so this ensures coverage
      // - these are loop limits (exclusive)
      const endX = Math.min(dimX, x + p + 2);
      const endY = Math.min(dimY, y + p + 2);
      const endZ = Math.min(dimZ, z + p + 2);

      // z/y/x with x, the texture-fast axis, innermost; each cell gets at most one contribution per vertex, so the
      // per-cell summation order is unchanged. `dz2 > p2` and `dy2 + dz2 > p2` only skip cells whose distance below
      // exceeds p as well (rounded addition of non-negative terms is monotone), and `d2` keeps the original
      // association so `d` is bit-identical; cells at `d2 > p2` but `d == p` would add `s = 0`, which changes nothing.
      for (let zi = begZ; zi < endZ; ++zi) {
        const dz = zi - vz;
        const dz2 = dz * dz;
        if (dz2 > p2) continue;
        const zo = zOff[zi];
        for (let yi = begY; yi < endY; ++yi) {
          const dy = yi - vy;
          const dy2 = dy * dy;
          if (dy2 + dz2 > p2) continue;
          const yo = zo + yi * width;
          for (let xi = begX; xi < endX; ++xi) {
            const dx = xi - vx;
            const d2 = dx * dx + dy2 + dz2;
            if (d2 > p2) continue;

            const s = p - Math.sqrt(d2);
            const c = yo + xi;
            if (itemSize === 4) {
              const index = c * 4;
              data[index] += c0 * s;
              data[index + 1] += c1 * s;
              data[index + 2] += c2 * s;
              data[index + 3] += c3 * s;
            } else {
              data[c] += c0 * s;
            }
            count[c] += s;
          }
        }
      }
    }
  }

  if (itemSize === 4) {
    for (let i = 0, il = count.length; i < il; ++i) {
      const is = i * 4;
      const c = count[i];
      grid[is] = Math.round(data[is] / c);
      grid[is + 1] = Math.round(data[is + 1] / c);
      grid[is + 2] = Math.round(data[is + 2] / c);
      grid[is + 3] = Math.round(data[is + 3] / c);
    }
  } else {
    for (let i = 0, il = count.length; i < il; ++i) {
      grid[i] = Math.round(data[i] / count[i]);
    }
  }

  const gridTexDim = Vec2.create(width, height);
  const gridTransform = Vec4.create(min[0], min[1], min[2], scaleFactor);
  const type = isInstanceType ? ('volumeInstance' as const) : ('volume' as const);

  if (graphics) {
    if (!texture) {
      const format = itemSize === 4 ? 'rgba' : 'r';
      texture = graphics.resources.texture('image-uint8', format, 'ubyte', 'linear');
    }
    texture.load(textureImage);

    return { kind: 'volume' as const, texture, gridTexDim, gridDim, gridTransform, type };
  } else {
    const interpolated = getTrilinearlyInterpolated({
      vertexCount,
      instanceCount,
      transformBuffer,
      positionBuffer,
      colorType: type,
      grid,
      gridDim,
      gridTexDim,
      gridTransform,
      vertexStride: 3,
      colorStride: itemSize,
      outputStride: itemSize,
    });

    return {
      kind: 'vertex' as const,
      texture: interpolated,
      texDim: Vec2.create(interpolated.width, interpolated.height),
      type: isInstanceType ? 'vertexInstance' : 'vertex',
    };
  }
}

/**
 * Offset of slice `z` in the 2d layout of a grid (`getVolumeTexture2dLayout`), in cells: cell (x, y, z) is at
 * `offsets[z] + y * width + x`
 */
function getSliceOffsets(gridDim: Vec3, width: number) {
  const [xn, yn, zn] = gridDim;
  const offsets = new Float64Array(zn);
  for (let z = 0; z < zn; ++z) {
    const column = Math.floor(((z * xn) % width) / xn);
    const row = Math.floor((z * xn) / width);
    offsets[z] = row * yn * width + column * xn;
  }
  return offsets;
}

//

interface ColorInterpolationInput {
  vertexCount: number;
  instanceCount: number;
  transformBuffer: Float32Array;
  positionBuffer: Float32Array;
  colorType: 'volumeInstance' | 'volume';
  grid: Uint8Array; // 2d layout
  gridTexDim: Vec2;
  gridDim: Vec3;
  gridTransform: Vec4;
  vertexStride: 3 | 4;
  colorStride: 1 | 3 | 4;
  outputStride: 1 | 3 | 4;
  itemOffset?: 0 | 1 | 2 | 3;
}

export function getTrilinearlyInterpolated(input: ColorInterpolationInput): TextureImage<Uint8Array> {
  const {
    vertexCount,
    positionBuffer,
    transformBuffer,
    grid,
    gridDim,
    gridTexDim,
    gridTransform,
    vertexStride,
    colorStride,
  } = input;

  const itemOffset = input.itemOffset || 0;
  const outputStride = input.outputStride;
  if (outputStride + itemOffset > colorStride) {
    throw new Error('outputStride + itemOffset must NOT be larger than colorStride');
  }

  const isInstanceType = input.colorType.endsWith('Instance');
  const instanceCount = isInstanceType ? input.instanceCount : 1;
  const image = createTextureImage(Math.max(1, instanceCount * vertexCount), outputStride, Uint8Array);
  const { array } = image;

  const [xn, yn, zn] = gridDim;
  const width = gridTexDim[0];
  const min = Vec3.fromArray(Vec3(), gridTransform, 0);
  const scaleFactor = gridTransform[3];

  const zOff = getSliceOffsets(gridDim, width);
  function getIndex(x: number, y: number, z: number) {
    // vertices outside the grid (not expected) keep their former, layout-derived index
    const zo =
      z >= 0 && z < zn ? zOff[z] : Math.floor(((z * xn) % width) / xn) * xn + Math.floor((z * xn) / width) * yn * width;
    return colorStride * (zo + y * width + x);
  }

  const v = Vec3();
  const v0 = Vec3();
  const v1 = Vec3();
  const vd = Vec3();

  for (let i = 0; i < instanceCount; ++i) {
    for (let j = 0; j < vertexCount; ++j) {
      Vec3.fromArray(v, positionBuffer, j * vertexStride);
      if (isInstanceType) Vec3.transformMat4Offset(v, v, transformBuffer, 0, 0, i * 16);
      Vec3.sub(v, v, min);
      Vec3.scale(v, v, scaleFactor);

      Vec3.floor(v0, v);
      Vec3.ceil(v1, v);

      Vec3.sub(vd, v, v0);
      Vec3.sub(v, v1, v0);
      Vec3.div(vd, vd, v);

      const [x0, y0, z0] = v0;
      const [x1, y1, z1] = v1;
      const [xd, yd, zd] = vd;

      const i000 = getIndex(x0, y0, z0) + itemOffset;
      const i100 = getIndex(x1, y0, z0) + itemOffset;
      const i001 = getIndex(x0, y0, z1) + itemOffset;
      const i101 = getIndex(x1, y0, z1) + itemOffset;
      const i010 = getIndex(x0, y1, z0) + itemOffset;
      const i110 = getIndex(x1, y1, z0) + itemOffset;
      const i011 = getIndex(x0, y1, z1) + itemOffset;
      const i111 = getIndex(x1, y1, z1) + itemOffset;

      const o = (i * vertexCount + j) * outputStride;

      for (let k = 0; k < outputStride; ++k) {
        const s000 = grid[i000 + k];
        const s100 = grid[i100 + k];
        const s001 = grid[i001 + k];
        const s101 = grid[i101 + k];
        const s010 = grid[i010 + k];
        const s110 = grid[i110 + k];
        const s011 = grid[i011 + k];
        const s111 = grid[i111 + k];

        const s00 = lerp(s000, s100, xd);
        const s01 = lerp(s001, s101, xd);
        const s10 = lerp(s010, s110, xd);
        const s11 = lerp(s011, s111, xd);

        const s0 = lerp(s00, s10, yd);
        const s1 = lerp(s01, s11, yd);

        array[o + k] = lerp(s0, s1, zd);
      }
    }
  }

  return image;
}

//

function isSupportedColorType(x: string): x is 'group' | 'groupInstance' {
  return x === 'group' || x === 'groupInstance';
}

export function applyMeshColorSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedColorType(values.colorType.ref.value)) return;

  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: values.aTransform.ref.value,
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tColor.ref.value,
      colorType: values.colorType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 4,
    },
    options,
    graphics,
    colorTexture,
  );

  if (smoothingData.kind === 'volume') {
    setColorType(values, smoothingData.type);
    ValueCell.update(values.tColorGrid, smoothingData.texture);
    ValueCell.update(values.uColorTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uColorGridDim, smoothingData.gridDim);
    ValueCell.update(values.uColorGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setColorType(values, smoothingData.type);
    ValueCell.update(values.tColor, smoothingData.texture);
    ValueCell.update(values.uColorTexDim, smoothingData.texDim);
  }
}

function isSupportedOverpaintType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshOverpaintSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedOverpaintType(values.overpaintType.ref.value)) return;

  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: values.aTransform.ref.value,
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tOverpaint.ref.value,
      colorType: values.overpaintType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 4,
    },
    options,
    graphics,
    colorTexture,
  );
  if (smoothingData.kind === 'volume') {
    setOverpaintType(values, smoothingData.type);
    ValueCell.update(values.tOverpaintGrid, smoothingData.texture);
    ValueCell.update(values.uOverpaintTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uOverpaintGridDim, smoothingData.gridDim);
    ValueCell.update(values.uOverpaintGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setOverpaintType(values, smoothingData.type);
    ValueCell.update(values.tOverpaint, smoothingData.texture);
    ValueCell.update(values.uOverpaintTexDim, smoothingData.texDim);
  }
}

function isSupportedTransparencyType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshTransparencySmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedTransparencyType(values.transparencyType.ref.value)) return;

  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: values.aTransform.ref.value,
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tTransparency.ref.value,
      colorType: values.transparencyType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 1,
    },
    options,
    graphics,
    colorTexture,
  );
  if (smoothingData.kind === 'volume') {
    setTransparencyType(values, smoothingData.type);
    ValueCell.update(values.tTransparencyGrid, smoothingData.texture);
    ValueCell.update(values.uTransparencyTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uTransparencyGridDim, smoothingData.gridDim);
    ValueCell.update(values.uTransparencyGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setTransparencyType(values, smoothingData.type);
    ValueCell.update(values.tTransparency, smoothingData.texture);
    ValueCell.update(values.uTransparencyTexDim, smoothingData.texDim);
  }
}

function isSupportedEmissiveType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshEmissiveSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedEmissiveType(values.emissiveType.ref.value)) return;

  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: values.aTransform.ref.value,
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tEmissive.ref.value,
      colorType: values.emissiveType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 1,
    },
    options,
    graphics,
    colorTexture,
  );
  if (smoothingData.kind === 'volume') {
    setEmissiveType(values, smoothingData.type);
    ValueCell.update(values.tEmissiveGrid, smoothingData.texture);
    ValueCell.update(values.uEmissiveTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uEmissiveGridDim, smoothingData.gridDim);
    ValueCell.update(values.uEmissiveGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setEmissiveType(values, smoothingData.type);
    ValueCell.update(values.tEmissive, smoothingData.texture);
    ValueCell.update(values.uEmissiveTexDim, smoothingData.texDim);
  }
}

function isSupportedSubstanceType(x: string): x is 'groupInstance' {
  return x === 'groupInstance';
}

export function applyMeshSubstanceSmoothing(
  values: MeshValues,
  options: ColorSmoothingOptions,
  graphics?: GraphicsContext,
  colorTexture?: Texture,
) {
  if (!isSupportedSubstanceType(values.substanceType.ref.value)) return;

  const smoothingData = calcMeshColorSmoothing(
    {
      vertexCount: values.uVertexCount.ref.value,
      instanceCount: values.uInstanceCount.ref.value,
      groupCount: values.uGroupCount.ref.value,
      transformBuffer: values.aTransform.ref.value,
      instanceBuffer: values.aInstance.ref.value,
      positionBuffer: values.aPosition.ref.value,
      groupBuffer: values.aGroup.ref.value,
      colorData: values.tSubstance.ref.value,
      colorType: values.substanceType.ref.value,
      boundingSphere: values.boundingSphere.ref.value,
      invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
      itemSize: 4,
    },
    options,
    graphics,
    colorTexture,
  );
  if (smoothingData.kind === 'volume') {
    setSubstanceType(values, smoothingData.type);
    ValueCell.update(values.tSubstanceGrid, smoothingData.texture);
    ValueCell.update(values.uSubstanceTexDim, smoothingData.gridTexDim);
    ValueCell.update(values.uSubstanceGridDim, smoothingData.gridDim);
    ValueCell.update(values.uSubstanceGridTransform, smoothingData.gridTransform);
  } else if (smoothingData.kind === 'vertex') {
    setSubstanceType(values, smoothingData.type);
    ValueCell.update(values.tSubstance, smoothingData.texture);
    ValueCell.update(values.uSubstanceTexDim, smoothingData.texDim);
  }
}
