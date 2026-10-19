/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util';
import { Box3D, Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
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
} from '../color-smoothing.js';
import { Mesh } from '../mesh.js';

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

/** smooths `layer` of a fresh mesh, on the CPU (per-vertex output) or for the GPU (grid output) */
function smooth(
  layer: Layer,
  type: 'group' | 'groupInstance',
  options: ColorSmoothingOptions,
  maxTextureSize?: number,
  loose = false,
) {
  const values = createValues(layer, type, loose);
  if (maxTextureSize === undefined) {
    apply[layer](values, options);
    const { type, texture, texDim } = layerValues(values, layer);
    return { type, texDim: [...texDim], hash: hashFnv32a(texture.array), array: texture.array };
  }
  const { graphics, loaded } = createGraphics(maxTextureSize);
  apply[layer](values, options, graphics);
  const { type: t, texDim, gridDim, gridTransform } = layerValues(values, layer);
  const image = loaded[0];
  return {
    type: t,
    texDim: [...texDim],
    gridDim: [...gridDim],
    gridTransform: [...gridTransform],
    hash: hashFnv32a(image.array),
    array: image.array,
  };
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

function summary(r: ReturnType<typeof smooth>) {
  const { array: _, ...rest } = r;
  return JSON.stringify(rest);
}

describe('mesh color smoothing', () => {
  it('pins the per-vertex output (no graphics context)', () => {
    const out: Record<string, string> = {};
    for (const [layer, type, loose] of Cases) {
      for (const o of Options) {
        out[`${caseName(layer, type, loose)} ${o.resolution}/${o.stride}`] = summary(
          smooth(layer, type, o, undefined, loose),
        );
      }
    }
    expect(out).toMatchInlineSnapshot(`
      {
        "color group 0.6/1": "{"type":"vertex","texDim":[56,54],"hash":2355055747}",
        "color group 1/3": "{"type":"vertex","texDim":[56,54],"hash":3257595857}",
        "color group 2.5/2": "{"type":"vertex","texDim":[56,54],"hash":360290621}",
        "color groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[96,94],"hash":1844755015}",
        "color groupInstance 1/3": "{"type":"vertexInstance","texDim":[96,94],"hash":1383235256}",
        "color groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[96,94],"hash":3399018012}",
        "emissive groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[95,95],"hash":1020822339}",
        "emissive groupInstance 1/3": "{"type":"vertexInstance","texDim":[95,95],"hash":1888119462}",
        "emissive groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[95,95],"hash":3316266964}",
        "loose color group 0.6/1": "{"type":"vertex","texDim":[56,54],"hash":360207778}",
        "loose color group 1/3": "{"type":"vertex","texDim":[56,54],"hash":925218032}",
        "loose color group 2.5/2": "{"type":"vertex","texDim":[56,54],"hash":1022845921}",
        "loose color groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[96,94],"hash":4259733340}",
        "loose color groupInstance 1/3": "{"type":"vertexInstance","texDim":[96,94],"hash":4201048997}",
        "loose color groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[96,94],"hash":4220855745}",
        "loose overpaint groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[96,94],"hash":1452890463}",
        "loose overpaint groupInstance 1/3": "{"type":"vertexInstance","texDim":[96,94],"hash":991663896}",
        "loose overpaint groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[96,94],"hash":1053528104}",
        "loose transparency groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[95,95],"hash":1318149860}",
        "loose transparency groupInstance 1/3": "{"type":"vertexInstance","texDim":[95,95],"hash":3832682863}",
        "loose transparency groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[95,95],"hash":330243369}",
        "overpaint groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[96,94],"hash":1833552715}",
        "overpaint groupInstance 1/3": "{"type":"vertexInstance","texDim":[96,94],"hash":1389030424}",
        "overpaint groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[96,94],"hash":378210256}",
        "substance groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[96,94],"hash":2282070990}",
        "substance groupInstance 1/3": "{"type":"vertexInstance","texDim":[96,94],"hash":3218325840}",
        "substance groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[96,94],"hash":1660208725}",
        "transparency groupInstance 0.6/1": "{"type":"vertexInstance","texDim":[95,95],"hash":4103678613}",
        "transparency groupInstance 1/3": "{"type":"vertexInstance","texDim":[95,95],"hash":2707318096}",
        "transparency groupInstance 2.5/2": "{"type":"vertexInstance","texDim":[95,95],"hash":887516233}",
      }
    `);
  });

  it('pins the grid output (graphics context), with and without a resolution coarsened to fit the texture', () => {
    const out: Record<string, string> = {};
    for (const [layer, type, loose] of Cases) {
      const name = caseName(layer, type, loose);
      for (const o of Options) {
        out[`${name} ${o.resolution}/${o.stride}`] = summary(smooth(layer, type, o, 16384, loose));
      }
      out[`${name} 0.6/1 max 128`] = summary(smooth(layer, type, Options[0], 128, loose));
    }
    expect(out).toMatchInlineSnapshot(`
      {
        "color group 0.6/1": "{"type":"volume","texDim":[240,155],"gridDim":[48,31,23],"gridTransform":[-10.598286247253418,-10.67664966583252,-1.1562506020069123,1.6666666666666667],"hash":2122876068}",
        "color group 0.6/1 max 128": "{"type":"volume","texDim":[124,84],"gridDim":[31,21,16],"gridTransform":[-10.964592247253417,-11.042955665832519,-1.5225566020069126,1.0348688717652583],"hash":1652453271}",
        "color group 1/3": "{"type":"volume","texDim":[120,84],"gridDim":[30,21,16],"gridTransform":[-10.998286247253418,-11.07664966583252,-1.5562506020069122,1],"hash":1640873575}",
        "color group 2.5/2": "{"type":"volume","texDim":[60,33],"gridDim":[15,11,9],"gridTransform":[-12.498286247253418,-12.57664966583252,-3.0562506020069122,0.4],"hash":598611831}",
        "color groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[497,370],"gridDim":[71,74,35],"gridTransform":[-10.598286247253418,-26.574967786924468,-1.1562506020069123,1.6666666666666667],"hash":800932756}",
        "color groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,124],"gridDim":[30,31,16],"gridTransform":[-11.554531723313419,-27.53121326298447,-2.1124960780669135,0.642572149049219],"hash":2192473603}",
        "color groupInstance 1/3": "{"type":"volumeInstance","texDim":[225,230],"gridDim":[45,46,23],"gridTransform":[-10.998286247253418,-26.974967786924466,-1.5562506020069122,1],"hash":2524297371}",
        "color groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,42],"gridDim":[21,21,12],"gridTransform":[-12.498286247253418,-28.474967786924466,-3.0562506020069122,0.4],"hash":2698519651}",
        "emissive groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[497,370],"gridDim":[71,74,35],"gridTransform":[-10.598286247253418,-26.574967786924468,-1.1562506020069123,1.6666666666666667],"hash":2949376451}",
        "emissive groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,124],"gridDim":[30,31,16],"gridTransform":[-11.554531723313419,-27.53121326298447,-2.1124960780669135,0.642572149049219],"hash":3763575608}",
        "emissive groupInstance 1/3": "{"type":"volumeInstance","texDim":[225,230],"gridDim":[45,46,23],"gridTransform":[-10.998286247253418,-26.974967786924466,-1.5562506020069122,1],"hash":1563449615}",
        "emissive groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,42],"gridDim":[21,21,12],"gridTransform":[-12.498286247253418,-28.474967786924466,-3.0562506020069122,0.4],"hash":2159287128}",
        "loose color group 0.6/1": "{"type":"volume","texDim":[240,160],"gridDim":[48,32,23],"gridTransform":[-10.8,-11.000000000000002,-1.4000000000000004,1.6666666666666667],"hash":624941653}",
        "loose color group 0.6/1 max 128": "{"type":"volume","texDim":[128,105],"gridDim":[32,21,17],"gridTransform":[-11.134775999999999,-11.303245999999996,-2.404327999999996,1.0348688717652583],"hash":395792937}",
        "loose color group 1/3": "{"type":"volume","texDim":[120,88],"gridDim":[30,22,16],"gridTransform":[-11,-12,-2,1],"hash":3288666885}",
        "loose color group 2.5/2": "{"type":"volume","texDim":[60,36],"gridDim":[15,12,9],"gridTransform":[-12.5,-15,-4,0.4],"hash":677397544}",
        "loose color groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[1008,320],"gridDim":[72,80,49],"gridTransform":[-10.8,-30.08770852684975,-4.4,1.6666666666666667],"hash":461035582}",
        "loose color groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,108],"gridDim":[24,27,17],"gridTransform":[-12.928637271364138,-32.544983069578024,-5.714549085456552,0.4827739662278129],"hash":849991221}",
        "loose color groupInstance 1/3": "{"type":"volumeInstance","texDim":[495,150],"gridDim":[45,50,31],"gridTransform":[-11,-30.687708526849747,-5,1],"hash":3998981198}",
        "loose color groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,72],"gridDim":[21,24,15],"gridTransform":[-12.5,-34.18770852684975,-6.5,0.4],"hash":746077568}",
        "loose overpaint groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[1008,320],"gridDim":[72,80,49],"gridTransform":[-10.8,-30.08770852684975,-4.4,1.6666666666666667],"hash":2048038893}",
        "loose overpaint groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,108],"gridDim":[24,27,17],"gridTransform":[-12.928637271364138,-32.544983069578024,-5.714549085456552,0.4827739662278129],"hash":3536861249}",
        "loose overpaint groupInstance 1/3": "{"type":"volumeInstance","texDim":[495,150],"gridDim":[45,50,31],"gridTransform":[-11,-30.687708526849747,-5,1],"hash":920044561}",
        "loose overpaint groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,72],"gridDim":[21,24,15],"gridTransform":[-12.5,-34.18770852684975,-6.5,0.4],"hash":493930730}",
        "loose transparency groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[1008,320],"gridDim":[72,80,49],"gridTransform":[-10.8,-30.08770852684975,-4.4,1.6666666666666667],"hash":3577676271}",
        "loose transparency groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,108],"gridDim":[24,27,17],"gridTransform":[-12.928637271364138,-32.544983069578024,-5.714549085456552,0.4827739662278129],"hash":1108770472}",
        "loose transparency groupInstance 1/3": "{"type":"volumeInstance","texDim":[495,150],"gridDim":[45,50,31],"gridTransform":[-11,-30.687708526849747,-5,1],"hash":2148475396}",
        "loose transparency groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,72],"gridDim":[21,24,15],"gridTransform":[-12.5,-34.18770852684975,-6.5,0.4],"hash":2615218260}",
        "overpaint groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[497,370],"gridDim":[71,74,35],"gridTransform":[-10.598286247253418,-26.574967786924468,-1.1562506020069123,1.6666666666666667],"hash":3747949446}",
        "overpaint groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,124],"gridDim":[30,31,16],"gridTransform":[-11.554531723313419,-27.53121326298447,-2.1124960780669135,0.642572149049219],"hash":3126137542}",
        "overpaint groupInstance 1/3": "{"type":"volumeInstance","texDim":[225,230],"gridDim":[45,46,23],"gridTransform":[-10.998286247253418,-26.974967786924466,-1.5562506020069122,1],"hash":504842673}",
        "overpaint groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,42],"gridDim":[21,21,12],"gridTransform":[-12.498286247253418,-28.474967786924466,-3.0562506020069122,0.4],"hash":59290004}",
        "substance groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[497,370],"gridDim":[71,74,35],"gridTransform":[-10.598286247253418,-26.574967786924468,-1.1562506020069123,1.6666666666666667],"hash":236527209}",
        "substance groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,124],"gridDim":[30,31,16],"gridTransform":[-11.554531723313419,-27.53121326298447,-2.1124960780669135,0.642572149049219],"hash":2046141351}",
        "substance groupInstance 1/3": "{"type":"volumeInstance","texDim":[225,230],"gridDim":[45,46,23],"gridTransform":[-10.998286247253418,-26.974967786924466,-1.5562506020069122,1],"hash":1028567617}",
        "substance groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,42],"gridDim":[21,21,12],"gridTransform":[-12.498286247253418,-28.474967786924466,-3.0562506020069122,0.4],"hash":3997557360}",
        "transparency groupInstance 0.6/1": "{"type":"volumeInstance","texDim":[497,370],"gridDim":[71,74,35],"gridTransform":[-10.598286247253418,-26.574967786924468,-1.1562506020069123,1.6666666666666667],"hash":1772353791}",
        "transparency groupInstance 0.6/1 max 128": "{"type":"volumeInstance","texDim":[120,124],"gridDim":[30,31,16],"gridTransform":[-11.554531723313419,-27.53121326298447,-2.1124960780669135,0.642572149049219],"hash":495575241}",
        "transparency groupInstance 1/3": "{"type":"volumeInstance","texDim":[225,230],"gridDim":[45,46,23],"gridTransform":[-10.998286247253418,-26.974967786924466,-1.5562506020069122,1],"hash":3598997994}",
        "transparency groupInstance 2.5/2": "{"type":"volumeInstance","texDim":[126,42],"gridDim":[21,21,12],"gridTransform":[-12.498286247253418,-28.474967786924466,-3.0562506020069122,0.4],"hash":3246908135}",
      }
    `);
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
