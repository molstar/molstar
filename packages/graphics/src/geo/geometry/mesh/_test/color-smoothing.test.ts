/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { ValueCell } from '@molstar/core/util';
import { Color } from '@molstar/core/util/color';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import {
  calculateTransformBoundingSphere,
  createTextureImage,
  type TextureImage,
} from '@molstar/graphics/gl/renderable/util';
import type { Texture } from '@molstar/graphics/gl/texture';
import { getVolumeTexture2dLayout } from '@molstar/graphics/repr/volume/util';
import { setColorType } from '../../color-data.js';
import { setEmissiveType } from '../../emissive-data.js';
import { setOverpaintType } from '../../overpaint-data.js';
import { setSubstanceType } from '../../substance-data.js';
import { createTransform } from '../../transform-data.js';
import { setTransparencyType } from '../../transparency-data.js';
import {
  applyMeshColorSmoothing,
  applyMeshEmissiveSmoothing,
  applyMeshOverpaintSmoothing,
  applyMeshSubstanceSmoothing,
  applyMeshTransparencySmoothing,
  type ColorSmoothingOptions,
  getTrilinearlyInterpolated,
} from '../color-smoothing.js';
import { Mesh } from '../mesh.js';
import { type ReferenceInput, referenceGrid, referenceInterpolate } from './color-smoothing-reference.js';

const VertexCount = 3000;
const GroupCount = 40;

function rng(seed: number) {
  return () => (seed = (seed * 16807) % 2147483647) / 2147483647;
}

/**
 * a lumpy ellipsoid shell off the origin; `loose` gives it a bounding sphere well beyond the vertices and without
 * extrema (surfaces get theirs from the structure), so its box is a cube much larger than the vertex box
 */
function createShell(loose: boolean) {
  const rnd = rng(7);
  const vertices = new Float32Array(VertexCount * 3);
  const groups = new Float32Array(VertexCount);
  for (let i = 0; i < VertexCount; ++i) {
    const u = rnd() * 2 - 1;
    const t = rnd() * Math.PI * 2;
    const s = Math.sqrt(1 - u * u);
    const r = 1 + 0.15 * Math.sin(5 * t) * u;
    vertices[i * 3] = 3 + 12 * r * s * Math.cos(t);
    vertices[i * 3 + 1] = -2 + 7 * r * s * Math.sin(t);
    vertices[i * 3 + 2] = 5 + 4 * r * u;
    groups[i] = Math.floor(((t / (2 * Math.PI)) * 8) | 0) * 5 + Math.min(4, Math.floor((u + 1) * 2.5));
  }
  const mesh = Mesh.create(vertices, new Uint32Array(0), new Float32Array(VertexCount * 3), groups, VertexCount, 0);
  if (loose) mesh.setBoundingSphere(Sphere3D.create(Vec3.create(4, -1, 5), 18));
  return mesh;
}

function createTransforms() {
  const ms = [
    Mat4.identity(),
    Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(20, 0, 0)),
      Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.create(0, 0, 1)),
    ),
    Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(0, -15, 8)),
      Mat4.fromRotation(Mat4(), (40 / 180) * Math.PI, Vec3.normalize(Vec3(), Vec3.create(1, 1, 0))),
    ),
  ];
  const array = new Float32Array(ms.length * 16);
  ms.forEach((m, i) => Mat4.toArray(m, array, i * 16));
  return createTransform(array, ms.length, undefined, 0, 0);
}

function fillImage(count: number, itemSize: 1 | 4, seed: number) {
  const image = createTextureImage(count, itemSize, Uint8Array);
  const rnd = rng(seed);
  for (let i = 0; i < count * itemSize; ++i) image.array[i] = Math.floor(rnd() * 256);
  return image;
}

type Layer = 'color' | 'overpaint' | 'transparency' | 'emissive' | 'substance';

function createValues(layer: Layer, type: 'group' | 'groupInstance', loose: boolean) {
  const mesh = createShell(loose);
  const transform = createTransforms();
  const values = Mesh.Utils.createValuesSimple(mesh, {}, Color(0x808080), 1, transform);
  ValueCell.update(values.uGroupCount, GroupCount);
  const instanceCount = values.uInstanceCount.ref.value;
  const count = type === 'group' ? GroupCount : GroupCount * instanceCount;
  switch (layer) {
    case 'color':
      ValueCell.update(values.tColor, fillImage(count, 4, 11));
      setColorType(values, type);
      break;
    case 'overpaint':
      ValueCell.update(values.tOverpaint, fillImage(count, 4, 13));
      setOverpaintType(values, type);
      break;
    case 'transparency':
      ValueCell.update(values.tTransparency, fillImage(count, 1, 17));
      setTransparencyType(values, type);
      break;
    case 'emissive':
      ValueCell.update(values.tEmissive, fillImage(count, 1, 19));
      setEmissiveType(values, type);
      break;
    case 'substance':
      ValueCell.update(values.tSubstance, fillImage(count, 4, 23));
      setSubstanceType(values, type);
      break;
  }
  return values;
}

/** just enough of a graphics context for the volume path: captures the image loaded into the grid texture */
function createGraphics(maxTextureSize: number) {
  const loaded: TextureImage<Uint8Array>[] = [];
  const graphics = {
    resources: {
      ownsTexture: () => true,
      limits: { maxTextureSize },
      texture: () => ({ load: (image: TextureImage<Uint8Array>) => loaded.push(image) }) as unknown as Texture,
    },
  } as unknown as GraphicsContext;
  return { graphics, loaded };
}

const apply = {
  color: applyMeshColorSmoothing,
  overpaint: applyMeshOverpaintSmoothing,
  transparency: applyMeshTransparencySmoothing,
  emissive: applyMeshEmissiveSmoothing,
  substance: applyMeshSubstanceSmoothing,
};

function layerValues(values: MeshValues, layer: Layer) {
  const prefix = layer === 'color' ? 'Color' : layer[0].toUpperCase() + layer.slice(1);
  const v = values as unknown as Record<string, ValueCell<any>>;
  return {
    type: v[`${layer}Type`].ref.value as string,
    texture: v[`t${prefix}`].ref.value as TextureImage<Uint8Array>,
    texDim: v[`u${prefix}TexDim`].ref.value as number[],
    gridDim: v[`u${prefix}GridDim`].ref.value as number[],
    gridTransform: v[`u${prefix}GridTransform`].ref.value as number[],
  };
}

/** the reference's input, read from `values` before smoothing replaces the per-vertex colors */
function referenceInput(values: MeshValues, layer: Layer, type: 'group' | 'groupInstance'): ReferenceInput {
  const isInstanceType = type === 'groupInstance';
  return {
    vertexCount: values.uVertexCount.ref.value,
    instanceCount: values.uInstanceCount.ref.value,
    groupCount: values.uGroupCount.ref.value,
    transformBuffer: values.aTransform.ref.value,
    positionBuffer: values.aPosition.ref.value,
    groupBuffer: values.aGroup.ref.value,
    colors: layerValues(values, layer).texture.array.slice(),
    isInstanceType,
    sphere: isInstanceType ? values.boundingSphere.ref.value : values.invariantBoundingSphere.ref.value,
    itemSize: layer === 'transparency' || layer === 'emissive' ? 1 : 4,
  };
}

/**
 * smooths `layer` of a fresh mesh, on the CPU (`array`: per-vertex output) or for the GPU (`array`: the grid), and
 * the reference's per-vertex output; `atVertices` is the output at the vertices either way
 */
function smooth(
  layer: Layer,
  type: 'group' | 'groupInstance',
  options: ColorSmoothingOptions,
  maxTextureSize?: number,
  loose = false,
) {
  const values = createValues(layer, type, loose);
  const input = referenceInput(values, layer, type);
  if (maxTextureSize === undefined) {
    apply[layer](values, options);
    const { type: t, texture } = layerValues(values, layer);
    const reference = referenceInterpolate(input, referenceGrid(input, options.resolution, options.stride));
    // the texture is padded beyond the vertices
    const atVertices = texture.array.subarray(0, reference.length);
    return { type: t, array: texture.array, atVertices, reference };
  }
  const { graphics, loaded } = createGraphics(maxTextureSize);
  apply[layer](values, options, graphics);
  const { type: t, texDim, gridDim, gridTransform } = layerValues(values, layer);
  const image = loaded[0];
  const scaleFactor = gridTransform[3];
  const atVertices = referenceInterpolate(input, {
    grid: image.array,
    gridDim,
    width: texDim[0],
    min: gridTransform,
    scaleFactor,
  });
  // the reference at the resolution actually used (coarsened to fit the texture)
  const referenceGridData = referenceGrid(input, 1 / scaleFactor, options.stride);
  return {
    type: t,
    array: image.array,
    image,
    texDim: [...texDim],
    gridDim: [...gridDim],
    resolution: 1 / scaleFactor,
    referenceGrid: referenceGridData,
    atVertices,
    reference: referenceInterpolate(input, referenceGridData),
  };
}

/**
 * at most 1 unit off, and in at most 1% of the values: floating-point reordering flips the odd rounding, a
 * systematic change (a bias, another kernel) moves many values
 */
function expectCloseTo(a: ArrayLike<number>, b: ArrayLike<number>, name: string) {
  expect(a.length, name).toBe(b.length);
  let max = 0;
  let differing = 0;
  for (let i = 0; i < a.length; ++i) {
    const d = Math.abs(a[i] - b[i]);
    if (d > 0) differing += 1;
    if (d > max) max = d;
  }
  expect(max, name).toBeLessThanOrEqual(1);
  expect(differing / a.length, name).toBeLessThanOrEqual(0.01);
}

/** not a degenerate output (all zero, constant): smoothed random group colors take many values */
function distinctValues(a: ArrayLike<number>) {
  return new Set(Array.from(a)).size;
}

const Cases: [Layer, 'group' | 'groupInstance', loose: boolean][] = [
  ['color', 'group', false],
  ['color', 'groupInstance', false],
  ['overpaint', 'groupInstance', false],
  ['transparency', 'groupInstance', false],
  ['emissive', 'groupInstance', false],
  ['substance', 'groupInstance', false],
  ['color', 'group', true],
  ['color', 'groupInstance', true],
  ['overpaint', 'groupInstance', true],
  ['transparency', 'groupInstance', true],
];

const caseName = (layer: Layer, type: string, loose: boolean) => `${loose ? 'loose ' : ''}${layer} ${type}`;

const Options: ColorSmoothingOptions[] = [
  { resolution: 0.6, stride: 1 },
  { resolution: 1, stride: 3 },
  { resolution: 2.5, stride: 2 },
];

describe('mesh color smoothing', () => {
  // the optimized kernel is bit-identical to the reference and the smaller grid keeps the reference's cell lattice;
  // `expectCloseTo` leaves room for floating-point reordering only
  it('matches the reference at the vertices (no graphics context)', () => {
    for (const [layer, type, loose] of Cases) {
      for (const o of Options) {
        const r = smooth(layer, type, o, undefined, loose);
        const name = `${caseName(layer, type, loose)} ${o.resolution}/${o.stride}`;
        expect(r.type, name).toBe(type === 'group' ? 'vertex' : 'vertexInstance');
        expectCloseTo(r.atVertices, r.reference, name);
        expect(distinctValues(r.atVertices), name).toBeGreaterThan(20);
      }
    }
  });

  it('gives a grid that matches the reference at the vertices, no larger than the bounding-sphere grid', () => {
    for (const [layer, type, loose] of Cases) {
      for (const o of Options) {
        const r = smooth(layer, type, o, 16384, loose);
        const name = `${caseName(layer, type, loose)} ${o.resolution}/${o.stride}`;
        expect(r.type, name).toBe(type === 'group' ? 'volume' : 'volumeInstance');
        expect(r.resolution, name).toBeCloseTo(o.resolution, 12);
        const { width, height } = getVolumeTexture2dLayout(Vec3.fromArray(Vec3(), r.gridDim, 0));
        expect(r.texDim, name).toEqual([width, height]);
        expect(r.image.array.length, name).toBe(
          width * height * (layer === 'transparency' || layer === 'emissive' ? 1 : 4),
        );
        for (let k = 0; k < 3; ++k) expect(r.gridDim[k], name).toBeLessThanOrEqual(r.referenceGrid.gridDim[k]);
        // a bounding sphere far beyond the vertices: the grid follows the vertices
        if (loose) expect(width * height, name).toBeLessThan(0.6 * r.referenceGrid.width * r.referenceGrid.height);
        expectCloseTo(r.atVertices, r.reference, name);
        expect(distinctValues(r.atVertices), name).toBeGreaterThan(20);
      }
    }
  });

  it('coarsens the resolution to fit the texture size', () => {
    for (const [layer, type, loose] of Cases) {
      const r = smooth(layer, type, Options[0], 128, loose);
      const name = caseName(layer, type, loose);
      expect(Math.max(...r.texDim), name).toBeLessThanOrEqual(128);
      expect(r.resolution, name).toBeGreaterThan(Options[0].resolution);
      expectCloseTo(r.atVertices, r.reference, name);
    }
  });

  it('gives the same output when run again, also after a larger grid', () => {
    const small = smooth('color', 'group', Options[2], 16384);
    const large = smooth('color', 'groupInstance', Options[0], 16384);
    expect(large.array.length).toBeGreaterThan(small.array.length);
    expect(smooth('color', 'group', Options[2], 16384).array).toEqual(small.array);
    expect(smooth('transparency', 'groupInstance', Options[1]).array).toEqual(
      smooth('transparency', 'groupInstance', Options[1]).array,
    );
  });

  it('spans the vertices of every instance with a margin, also on a grid much smaller than the bounding-sphere box', () => {
    // a third instance with a reflection and a non-uniform scale
    const skewed = Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(-30, 4, 2)),
      Mat4.fromScaling(Mat4(), Vec3.create(-1, 1.7, 0.6)),
    );
    for (const [layer, type, loose] of Cases) {
      for (const skew of type === 'groupInstance' ? [false, true] : [false]) {
        for (const o of [...Options, { resolution: 0.25, stride: 1 }]) {
          const values = createValues(layer, type, loose);
          const transforms = values.aTransform.ref.value;
          if (skew) {
            Mat4.toArray(skewed, transforms, 2 * 16);
            const bs = calculateTransformBoundingSphere(values.invariantBoundingSphere.ref.value, transforms, 3, 0);
            ValueCell.update(values.boundingSphere, bs);
          }
          const { graphics } = createGraphics(16384);
          apply[layer](values, o, graphics);
          const { gridDim, gridTransform } = layerValues(values, layer);
          const positions = values.aPosition.ref.value;
          const instanceCount = type === 'groupInstance' ? values.uInstanceCount.ref.value : 1;
          const sphere = type === 'groupInstance' ? values.boundingSphere : values.invariantBoundingSphere;
          const sphereBox = Box3D.fromSphere3D(Box3D(), sphere.ref.value);
          const p = Vec3();
          let lo = Infinity;
          let hi = -Infinity;
          for (let i = 0; i < instanceCount; ++i) {
            for (let j = 0; j < VertexCount; ++j) {
              Vec3.fromArray(p, positions, j * 3);
              if (type === 'groupInstance') Vec3.transformMat4Offset(p, p, transforms, 0, 0, i * 16);
              for (let k = 0; k < 3; ++k) {
                const g = (p[k] - gridTransform[k]) * gridTransform[3];
                // a vertex outside the bounding-sphere box (the box of transformed extrema can miss some) loses that
                // much margin, as on the grid of the bounding-sphere box
                lo = Math.min(lo, g + Math.max(0, (sphereBox.min[k] - p[k]) * gridTransform[3]));
                hi = Math.max(hi, g - gridDim[k] - Math.max(0, (p[k] - sphereBox.max[k]) * gridTransform[3]));
              }
            }
          }
          // the box is expanded by 1 + resolution (1 + 1 / resolution cells) and has 2 extra cells at the end
          const padCells = 1 + 1 / o.resolution;
          expect(lo).toBeGreaterThanOrEqual(padCells - 1e-6);
          expect(hi).toBeLessThanOrEqual(-(2 + padCells) + 1e-6);
        }
      }
    }
  });

  it('smooths a mesh without vertices', () => {
    const values = createValues('color', 'groupInstance', false);
    ValueCell.update(values.uVertexCount, 0);
    const { graphics, loaded } = createGraphics(16384);
    applyMeshColorSmoothing(values, Options[1], graphics);
    expect(loaded[0].array.every((x) => x === 0)).toBe(true);
    applyMeshColorSmoothing(createValues('color', 'group', false), Options[1]);
  });
});

describe('getTrilinearlyInterpolated', () => {
  it('interpolates at vertices on grid planes (perf-audit-2 08/B3)', () => {
    // a uniform 3x3x3 grid of value 200, one unit per cell, its origin at 0
    const gridDim = Vec3.create(3, 3, 3);
    const { width, height } = getVolumeTexture2dLayout(gridDim);
    const grid = new Uint8Array(width * height).fill(200);
    // on a corner, on an edge, on a face, inside a cell
    const positionBuffer = new Float32Array([1, 1, 1, 1, 1, 0.5, 1, 0.5, 0.5, 0.5, 0.5, 0.5]);
    const image = getTrilinearlyInterpolated({
      vertexCount: 4,
      instanceCount: 1,
      transformBuffer: new Float32Array(Mat4.identity()),
      positionBuffer,
      colorType: 'volume',
      grid,
      gridTexDim: Vec2.create(width, height),
      gridDim,
      gridTransform: Vec4.create(0, 0, 0, 1),
      vertexStride: 3,
      colorStride: 1,
      outputStride: 1,
    });
    expect(Array.from(image.array.subarray(0, 4))).toEqual([200, 200, 200, 200]);
  });
});
