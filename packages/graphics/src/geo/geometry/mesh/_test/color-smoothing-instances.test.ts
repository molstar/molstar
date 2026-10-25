/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import { hashFnv32a } from '@molstar/core/data/util/hash-functions';
import { Sphere3D } from '@molstar/core/math/geometry';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { RuntimeContext } from '@molstar/core/task';
import { ValueCell } from '@molstar/core/util';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { CIF } from '@molstar/io/reader/cif';
import { trajectoryFromMmCIF } from '@molstar/model/formats/structure/mmcif';
import { Structure } from '@molstar/model/model/structure';
import { StructureSymmetry } from '@molstar/model/model/structure/structure/symmetry';
import type { CustomProperty } from '@molstar/model/props/common/custom-property';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import type { MeshValues } from '@molstar/graphics/gl/renderable/mesh';
import { createDataArray, type DataArray, type TextureImage } from '@molstar/graphics/gl/renderable/util';
import type { Texture } from '@molstar/graphics/gl/texture';
import type { Representation } from '@molstar/graphics/repr/representation';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { ColorTheme } from '@molstar/graphics/theme/color';
import { BuiltInColorThemes } from '@molstar/graphics/theme/color/catalog';
import { SizeTheme } from '@molstar/graphics/theme/size';
import { BuiltInSizeThemes } from '@molstar/graphics/theme/size/catalog';
import { Theme } from '@molstar/graphics/theme/theme';
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
  calcMeshColorSmoothing,
} from '../color-smoothing.js';
import { Mesh } from '../mesh.js';

// Layer smoothing on an instanced surface whose layer reaches only some instances (perf-audit 13/P1): the gaussian
// surface of the 1cbs protein over two unit cells (8 instances, neighbours within the smoothing reach). The snapshot
// pins, per instance, what the CPU path writes and what the GPU path's grid gives at the vertices, as computed when
// every instance was splatted over the grid of all instances.

const registries = (() => {
  const colorThemeRegistry = ColorTheme.createRegistry();
  for (const p of Object.values(BuiltInColorThemes)) colorThemeRegistry.add(p as ColorTheme.Provider<any, any, any>);
  const sizeThemeRegistry = SizeTheme.createRegistry();
  for (const p of Object.values(BuiltInSizeThemes)) sizeThemeRegistry.add(p as SizeTheme.Provider<any, any>);
  return { colorThemeRegistry, sizeThemeRegistry };
})();
const propCtx: CustomProperty.Context = {
  runtime: RuntimeContext.Synchronous,
  assetManager: new DefaultAssetManager(),
};

async function loadCells() {
  const data = fs.readFileSync(
    path.resolve(import.meta.dirname, '../../../../../../../data/examples/1cbs_updated.cif'),
    'utf8',
  );
  const parsed = await CIF.parse(data).run();
  if (parsed.isError) throw new Error(parsed.message);
  const trajectory = await trajectoryFromMmCIF(parsed.result.blocks[0], parsed.result).run();
  const model = Structure.ofModel(trajectory.representative);
  return StructureSymmetry.buildSymmetryRange(model, Vec3.create(0, 0, 0), Vec3.create(1, 0, 0)).run();
}

/** the geometry of the largest gaussian-surface visual: positions, groups, transforms and bounding sphere */
async function surfaceGeometry(structure: Structure) {
  const provider = GaussianSurfaceRepresentationProvider;
  const themeProps = {
    colorTheme: { name: provider.defaultColorTheme.name, params: provider.defaultColorTheme.props ?? {} },
    sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
  };
  await Theme.ensureDependencies(propCtx, registries, { structure }, themeProps);
  const repr = provider.factory({ ...registries }, provider.getParams) as Representation<Structure>;
  repr.setTheme(Theme.create(registries, { structure }, themeProps));
  const params = PD.getDefaultValues(provider.getParams(registries, structure));
  await repr
    .createOrUpdate({ ...params, renderMerged: false, smoothColors: { name: 'off', params: {} } }, structure)
    .run();
  let values: MeshValues | undefined;
  for (const o of repr.renderObjects) {
    const v = o.values as MeshValues;
    if (!values || v.uVertexCount.ref.value > values.uVertexCount.ref.value) values = v;
  }
  if (!values) throw new Error('expected a surface');
  const geometry = {
    vertexCount: values.uVertexCount.ref.value,
    groupCount: values.uGroupCount.ref.value,
    instanceCount: values.uInstanceCount.ref.value,
    positions: values.aPosition.ref.value.slice(0, values.uVertexCount.ref.value * 3),
    groups: values.aGroup.ref.value.slice(0, values.uVertexCount.ref.value),
    transforms: values.aTransform.ref.value.slice(0, values.uInstanceCount.ref.value * 16),
    invariantBoundingSphere: Sphere3D.clone(values.invariantBoundingSphere.ref.value),
  };
  repr.destroy();
  return geometry;
}

type Geometry = Awaited<ReturnType<typeof surfaceGeometry>> & {
  /** `aInstance` (and `aTransform` in its order), as reordered by the instance grid */
  permutation?: number[];
  /** the layer type, `groupInstance` when not given; `group` for color only */
  colorType?: 'group' | 'groupInstance';
  /** scales the bounding sphere: below 1 some vertices fall outside its box, and the grid */
  sphereScale?: number;
};
type Layer = 'color' | 'overpaint' | 'transparency' | 'emissive' | 'substance';

const ItemSize: Record<Layer, 1 | 4> = { color: 4, overpaint: 4, transparency: 1, emissive: 1, substance: 4 };

function rng(seed: number) {
  return () => (seed = (seed * 16807) % 2147483647) / 2147483647;
}

/** a `groupInstance` layer, non-zero on the groups `reach(instance, group)` accepts, random values there */
function layerImage(geometry: Geometry, itemSize: 1 | 4, reach: (instance: number, group: number) => boolean) {
  const { groupCount } = geometry;
  const instanceCount = geometry.colorType === 'group' ? 1 : geometry.instanceCount;
  const image = createDataArray(groupCount * instanceCount, itemSize, Uint8Array);
  const rnd = rng(31);
  for (let i = 0; i < instanceCount; ++i) {
    for (let g = 0; g < groupCount; ++g) {
      const o = (i * groupCount + g) * itemSize;
      for (let k = 0; k < itemSize; ++k) {
        const x = 1 + Math.floor(rnd() * 255);
        if (reach(i, g)) image.array[o + k] = x;
      }
    }
  }
  return image;
}

function createValues(geometry: Geometry, layer: Layer, image: DataArray<Uint8Array>) {
  const { positions, groups, vertexCount, transforms, instanceCount } = geometry;
  const mesh = Mesh.create(
    positions.slice(),
    new Uint32Array(0),
    new Float32Array(vertexCount * 3),
    groups.slice(),
    vertexCount,
    0,
  );
  const sphere = geometry.invariantBoundingSphere;
  // a scaled sphere has no extrema, so its box is the cube of the scaled radius
  const { sphereScale } = geometry;
  mesh.setBoundingSphere(sphereScale ? Sphere3D.create(sphere.center, sphere.radius * sphereScale) : sphere);
  const transform = createTransform(transforms.slice(), instanceCount, undefined, 0, 0);
  const values = Mesh.Utils.createValuesSimple(mesh, {}, Color(0x808080), 1, transform);
  ValueCell.update(values.uGroupCount, geometry.groupCount);
  const { permutation } = geometry;
  if (permutation) {
    const aInstance = values.aInstance.ref.value;
    const aTransform = values.aTransform.ref.value;
    permutation.forEach((p, k) => {
      aInstance[k] = p;
      aTransform.set(transforms.subarray(p * 16, (p + 1) * 16), k * 16);
    });
  }
  const type = geometry.colorType ?? 'groupInstance';
  switch (layer) {
    case 'color':
      ValueCell.update(values.tColor, image);
      setColorType(values, type);
      break;
    case 'overpaint':
      ValueCell.update(values.tOverpaint, image);
      setOverpaintType(values, 'groupInstance');
      break;
    case 'transparency':
      ValueCell.update(values.tTransparency, image);
      setTransparencyType(values, 'groupInstance');
      break;
    case 'emissive':
      ValueCell.update(values.tEmissive, image);
      setEmissiveType(values, 'groupInstance');
      break;
    case 'substance':
      ValueCell.update(values.tSubstance, image);
      setSubstanceType(values, 'groupInstance');
      break;
  }
  return values;
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
    texture: v[`t${prefix}`].ref.value as DataArray<Uint8Array>,
    texDim: v[`u${prefix}GridTexDim`].ref.value as number[],
    gridDim: v[`u${prefix}GridDim`].ref.value as number[],
    gridTransform: v[`u${prefix}GridTransform`].ref.value as number[],
  };
}

/** just enough of a graphics context for the volume path: captures the image loaded into the grid texture */
function createGraphics() {
  const loaded: TextureImage<Uint8Array>[] = [];
  const graphics = {
    resources: {
      ownsTexture: () => true,
      limits: { maxTextureSize: 16384 },
      texture: () => ({ load: (image: TextureImage<Uint8Array>) => loaded.push(image) }) as unknown as Texture,
    },
  } as unknown as GraphicsContext;
  return { graphics, loaded };
}

function hash(array: Uint8Array) {
  return hashFnv32a(array as unknown as number[]).toString(16);
}

function isZero(array: Uint8Array) {
  for (let i = 0; i < array.length; ++i) if (array[i]) return false;
  return true;
}

/** per instance (output row): `0` when all clear, else the hash of its values */
function perInstance(array: Uint8Array, geometry: Geometry, itemSize: number) {
  const n = geometry.vertexCount * itemSize;
  const out: string[] = [];
  const rows = geometry.colorType === 'group' ? 1 : geometry.instanceCount;
  for (let i = 0; i < rows; ++i) {
    const row = array.subarray(i * n, (i + 1) * n);
    out.push(isZero(row) ? '0' : hash(row));
  }
  return out.join(' ');
}

/** the CPU output (per vertex, per instance) and the GPU path's grid */
function smoothed(geometry: Geometry, layer: Layer, image: DataArray<Uint8Array>, options: ColorSmoothingOptions) {
  const itemSize = ItemSize[layer];
  const cpuValues = createValues(geometry, layer, { ...image, array: image.array.slice() });
  apply[layer](cpuValues, options);
  const cpu = layerValues(cpuValues, layer);

  const gpuValues = createValues(geometry, layer, { ...image, array: image.array.slice() });
  const { graphics, loaded } = createGraphics();
  apply[layer](gpuValues, options, graphics);
  const gpu = layerValues(gpuValues, layer);
  const grid = loaded[0].array;
  return {
    cpu: `${cpu.type} ${perInstance(cpu.texture.array, geometry, itemSize)}`,
    gpu: `${gpu.type} ${[...gpu.gridDim].join('x')} ${[...gpu.texDim].join('x')} ${[...gpu.gridTransform].map((x) => x.toFixed(4)).join(',')} grid:${isZero(grid) ? '0' : hash(grid)}`,
    cpuZero: isZero(cpu.texture.array),
    gridZero: isZero(grid),
  };
}

/** two pairs of neighbours within the smoothing reach (0 and 1; 2 and the rotated 3), then 4 and 5 farther out */
function spreadTransforms() {
  const ms = [
    Mat4.identity(),
    Mat4.fromTranslation(Mat4(), Vec3.create(30, 0, 0)),
    Mat4.fromTranslation(Mat4(), Vec3.create(0, 90, 0)),
    Mat4.mul(
      Mat4(),
      Mat4.fromTranslation(Mat4(), Vec3.create(25, 100, 5)),
      Mat4.fromRotation(Mat4(), Math.PI / 2, Vec3.create(0, 0, 1)),
    ),
    Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -90)),
    Mat4.fromTranslation(Mat4(), Vec3.create(60, 50, -45)),
  ];
  const transforms = new Float32Array(ms.length * 16);
  ms.forEach((m, i) => Mat4.toArray(m, transforms, i * 16));
  return transforms;
}

/** the surface, computed once for the file: the cases copy what they change */
const CellsGeometry = loadCells().then(surfaceGeometry);

/** the first layer's sweep is the quick lane's representative, the others are medium (AGENTS.md, Tests) */
const sweepLane = (layer: Layer) => (layer === 'overpaint' ? {} : { tags: ['medium'] });

const Options: ColorSmoothingOptions[] = [
  { resolution: 1.2, stride: 3 },
  { resolution: 0.8, stride: 2 },
];

describe('mesh layer smoothing on instances', async () => {
  const geometry = await CellsGeometry;
  const { groupCount, instanceCount } = geometry;

  const Reach: [string, (instance: number, group: number) => boolean, ColorSmoothingOptions[]][] = [
    ['one instance', (i) => i === 3, Options],
    ['first instance, some groups', (i, g) => i === 0 && g >= groupCount / 4 && g < groupCount / 2, Options],
    ['two instances', (i) => i === 1 || i === 6, [Options[0]]],
    ['every instance', () => true, [Options[0]]],
    ['none', () => false, [Options[0]]],
  ];

  it('has an instanced fixture', () => {
    expect(instanceCount).toBe(8);
    expect(geometry.vertexCount).toBeGreaterThan(1000);
  });

  for (const layer of ['overpaint', 'transparency', 'emissive', 'substance', 'color'] as Layer[]) {
    it(`smooths ${layer} reaching some instances as when every instance is splatted`, sweepLane(layer), () => {
      const result: Record<string, ReturnType<typeof smoothed>> = {};
      for (const [name, reach, options] of Reach) {
        for (const o of options) {
          const r = smoothed(geometry, layer, layerImage(geometry, ItemSize[layer], reach), o);
          result[`${name} ${o.resolution}/${o.stride}`] = r;
          if (name === 'none') {
            expect(r.cpuZero, name).toBe(true);
            expect(r.gridZero, name).toBe(true);
          } else {
            expect(r.cpuZero, name).toBe(false);
            expect(r.gridZero, name).toBe(false);
          }
        }
      }
      expect(result).toMatchSnapshot();
    });
  }
});

describe('mesh layer smoothing on spread instances', async () => {
  const cells = await CellsGeometry;
  const transforms = spreadTransforms();
  const geometry: Geometry = { ...cells, instanceCount: transforms.length / 16, transforms };
  const { groupCount } = geometry;

  const Reach: [string, (instance: number, group: number) => boolean, ColorSmoothingOptions[]][] = [
    ['one instance with a neighbour', (i) => i === 0, Options],
    ['the rotated instance', (i) => i === 3, [Options[0]]],
    ['a lone instance, some groups', (i, g) => i === 4 && g >= groupCount / 4 && g < groupCount / 2, Options],
    ['two lone instances', (i) => i === 4 || i === 5, [Options[0]]],
    ['none', () => false, [Options[0]]],
  ];

  for (const layer of ['overpaint', 'transparency', 'emissive', 'substance'] as Layer[]) {
    it(`smooths ${layer} reaching some instances as when every instance is splatted`, sweepLane(layer), () => {
      const result: Record<string, ReturnType<typeof smoothed>> = {};
      for (const [name, reach, options] of Reach) {
        for (const o of options) {
          const r = smoothed(geometry, layer, layerImage(geometry, ItemSize[layer], reach), o);
          result[`${name} ${o.resolution}/${o.stride}`] = r;
          expect(r.cpuZero, name).toBe(name === 'none');
          expect(r.gridZero, name).toBe(name === 'none');
        }
      }
      expect(result).toMatchSnapshot();
    });
  }
});

describe('mesh layer smoothing on reordered and clipped instances, odd groups', async () => {
  const cells = await CellsGeometry;
  const transforms = spreadTransforms();
  const spread: Geometry = { ...cells, instanceCount: transforms.length / 16, transforms };
  const outOfRange = spread.groups.slice();
  outOfRange[10] = spread.groupCount;
  outOfRange[20] = spread.groupCount + 3;
  const quarter = (g: number) => g >= spread.groupCount / 4 && g < spread.groupCount / 2;

  const Cases: [string, Layer, Geometry, (instance: number, group: number) => boolean][] = [
    // the GPU path reads the rows of `aInstance`, the CPU path the serial rows
    ['reordered, a lone instance', 'overpaint', { ...spread, permutation: [3, 5, 0, 4, 1, 2] }, (i) => i === 4],
    ['reordered, a pair', 'transparency', { ...spread, permutation: [3, 5, 0, 4, 1, 2] }, (i) => i === 0],
    ['clipped sphere, a pair', 'overpaint', { ...spread, sphereScale: 0.5 }, (i) => i === 0],
    ['clipped sphere, a lone instance', 'emissive', { ...spread, sphereScale: 0.5 }, (i) => i === 4],
    ['groups out of range', 'overpaint', { ...spread, groups: outOfRange }, (i) => i === 4],
    ['group type', 'color', { ...spread, colorType: 'group' }, (_, g) => quarter(g)],
    ['group type, nothing', 'color', { ...spread, colorType: 'group' }, () => false],
  ];

  it('smooths as when every instance is splatted', () => {
    const result: Record<string, ReturnType<typeof smoothed>> = {};
    for (const [name, layer, geometry, reach] of Cases) {
      result[name] = smoothed(geometry, layer, layerImage(geometry, ItemSize[layer], reach), Options[0]);
    }
    expect(result).toMatchSnapshot();
  });
});

describe('mesh layer smoothing splats the reached instances and their neighbours', async () => {
  const cells = await CellsGeometry;
  const transforms = spreadTransforms();
  const spread: Geometry = { ...cells, instanceCount: transforms.length / 16, transforms };

  /** the instances splatted for an overpaint reaching `instances`, without and with a graphics context */
  function splatCounts(geometry: Geometry, instances: number[]) {
    return [false, true].map((withGraphics) => {
      const values = createValues(
        geometry,
        'overpaint',
        layerImage(geometry, 4, (i) => instances.includes(i)),
      );
      const r = calcMeshColorSmoothing(
        {
          vertexCount: values.uVertexCount.ref.value,
          instanceCount: values.uInstanceCount.ref.value,
          groupCount: values.uGroupCount.ref.value,
          transformBuffer: values.aTransform.ref.value,
          instanceBuffer: values.aInstance.ref.value,
          positionBuffer: values.aPosition.ref.value,
          groupBuffer: values.aGroup.ref.value,
          colorData: values.tOverpaint.ref.value,
          colorType: 'groupInstance',
          boundingSphere: values.boundingSphere.ref.value,
          invariantBoundingSphere: values.invariantBoundingSphere.ref.value,
          itemSize: 4,
        },
        Options[0],
        withGraphics ? createGraphics().graphics : undefined,
      );
      return r.splatCount;
    });
  }

  it('skips instances out of reach', () => {
    expect(splatCounts(spread, [])).toEqual([0, 0]);
    expect(splatCounts(spread, [4])).toEqual([1, 1]);
    expect(splatCounts(spread, [0])).toEqual([2, 2]);
    // two reached instances far apart: their neighbours, not every instance between them
    expect(splatCounts(spread, [0, 4])).toEqual([3, 3]);
    expect(splatCounts(spread, [0, 1, 2, 3, 4, 5])).toEqual([6, 6]);
    // with a graphics context, rows are those of `aInstance`: instance 4 is drawn at index 3
    expect(splatCounts({ ...spread, permutation: [3, 5, 0, 4, 1, 2] }, [4])).toEqual([2, 1]);
  });
});
