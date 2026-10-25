/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type GraphicsContext } from '@molstar/graphics/gl/context';
import { Box3D } from '@molstar/core/math/geometry';
import { Mat4, Tensor, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import { ValueCell } from '@molstar/core/util';
import { createGrid3dComputeRenderable } from '@molstar/graphics/gl/compute/grid3d';
import { TextureSpec, UniformSpec } from '@molstar/graphics/gl/renderable/schema';
import { CustomProperties } from '@molstar/model/model/custom-property';
import { Volume } from '@molstar/model/model/volume';
import { Representation, RepresentationProvider } from '@molstar/graphics/repr/representation';
import { getThemeRegistries } from './structure.js';
import { Theme } from '@molstar/graphics/theme/theme';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';

/** Gaussian blobs [x, y, z, sigma] on an n^3 grid with unit spacing */
function getBlobs(n: number) {
  const c = n / 2;
  return [
    [c - n / 6, c, c, n / 5],
    [c + n / 6, c, c - n / 8, n / 6],
    [c, c + n / 5, c + n / 8, n / 8],
  ];
}

/** `volume` with negated values */
export function negatedVolume(volume: Volume): Volume {
  const { space, data } = volume.grid.cells;
  return createVolume(
    `negated ${volume.label}`,
    space,
    Float32Array.from(data as ArrayLike<number>, (v) => -v),
  );
}

/** Sum of three Gaussian blobs on an n^3 grid with unit spacing */
export function blobVolume(n: number): Volume {
  const space = Tensor.Space([n, n, n], [2, 1, 0], Float32Array);
  const data = new Float32Array(n * n * n);
  const blobs = getBlobs(n);
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) {
        let v = 0;
        for (const [x, y, z, s] of blobs) {
          const d2 = (i - x) ** 2 + (j - y) ** 2 + (k - z) ** 2;
          v += Math.exp(-d2 / (2 * s * s));
        }
        data[space.dataOffset(i, j, k)] = v;
      }
    }
  }
  return createVolume('blobs', space, data);
}

const BlobSum = createGrid3dComputeRenderable({
  schema: {
    tBlobs: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
    uBlobCount: UniformSpec('i'),
  },
  loopBounds: ['uBlobCount'],
  mainCode: `
    float v = 0.0;
    for (int b = 0; b < uBlobCount; ++b) {
      vec4 blob = texture2D(tBlobs, vec2((float(b) + 0.5) / float(uBlobCount), 0.5));
      vec3 d = xyz - blob.xyz;
      v += exp(-dot(d, d) / (2.0 * blob.w * blob.w));
    }`,
  returnCode: 'v',
  wgsl: {
    mainCode: /* wgsl */ `
    var v = 0.0;
    for (var b = 0; b < uParams.uBlobCount; b++) {
      let blob = textureLoad(tBlobs, vec2<i32>(b, 0), 0);
      let d = xyz - blob.xyz;
      v += exp(-dot(d, d) / (2.0 * blob.w * blob.w));
    }`,
    returnCode: 'v',
  },
  values: (blobs: number[][]) => ({
    tBlobs: { width: blobs.length, height: 1, array: new Float32Array(blobs.flat()) },
    uBlobCount: blobs.length,
  }),
});

const BlobSumCumulative = createGrid3dComputeRenderable({
  schema: { uBlob: UniformSpec('v4') },
  mainCode: 'vec3 d = xyz - uBlob.xyz;',
  returnCode: 'current + exp(-dot(d, d) / (2.0 * uBlob.w * uBlob.w))',
  wgsl: {
    mainCode: 'let d = xyz - uParams.uBlob.xyz;',
    returnCode: 'current + exp(-dot(d, d) / (2.0 * uParams.uBlob.w * uParams.uBlob.w))',
  },
  values: () => ({ uBlob: Vec4() }),
  cumulative: {
    states: (blobs: number[][]) => blobs,
    update: (_, blob, values) => ValueCell.update(values.uBlob, Vec4.fromArray(Vec4(), blob, 0)),
  },
});

/** Same as `blobVolume`, computed with `createGrid3dComputeRenderable` (in one pass or cumulative per blob) */
export async function gpuBlobVolume(graphics: GraphicsContext, n: number, cumulative: boolean): Promise<Volume> {
  const grid = { box: Box3D.create(Vec3(), Vec3.create(n - 1, n - 1, n - 1)), dimensions: Vec3.create(n, n, n) };
  const compute = cumulative ? BlobSumCumulative : BlobSum;
  const data = await compute(SyncRuntimeContext, graphics, grid, getBlobs(n));
  // grid3d axis order is fast to slow z, y, x
  return createVolume('gpu-blobs', Tensor.Space([n, n, n], [0, 1, 2], Float32Array), data);
}

/**
 * Blobs as a segmentation: each cell labeled by its strongest blob where that exceeds 0.5 (else 0), a segment per label
 * with the bounds of its cells (grid coordinates, inclusive), as segmentation volumes have them
 */
export function segmentedBlobVolume(n: number): Volume {
  const space = Tensor.Space([n, n, n], [2, 1, 0], Float32Array);
  const data = new Float32Array(n * n * n);
  const blobs = getBlobs(n);
  const bounds: { [k: Volume.SegmentIndex]: Box3D } = {};
  const segments = new Map<Volume.SegmentIndex, Set<number>>();
  const sets = new Map<number, Set<Volume.SegmentIndex>>();
  blobs.forEach((_, b) => {
    const segment = (b + 1) as Volume.SegmentIndex;
    bounds[segment] = Box3D.create(Vec3.create(n - 1, n - 1, n - 1), Vec3.create(-1, -1, -1));
    segments.set(segment, new Set([b + 1]));
    sets.set(b + 1, new Set([segment]));
  });
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) {
        let label = 0;
        let strongest = 0.5;
        blobs.forEach(([x, y, z, s], b) => {
          const v = Math.exp(-((i - x) ** 2 + (j - y) ** 2 + (k - z) ** 2) / (2 * s * s));
          if (v > strongest) {
            strongest = v;
            label = b + 1;
          }
        });
        data[space.dataOffset(i, j, k)] = label;
        if (label === 0) continue;
        const { min, max } = bounds[label as Volume.SegmentIndex];
        Vec3.min(min, min, Vec3.create(i, j, k));
        Vec3.max(max, max, Vec3.create(i, j, k));
      }
    }
  }
  const volume = createVolume('segmented blobs', space, data);
  Volume.Segmentation.set(volume, { segments, sets, bounds, labels: {} });
  Volume.PickingGranularity.set(volume, 'object');
  return volume;
}

export function createVolume(label: string, space: Tensor.Space, data: Float32Array): Volume {
  let min = Infinity,
    max = -Infinity,
    sum = 0;
  for (let i = 0; i < data.length; ++i) {
    const v = data[i];
    if (v < min) min = v;
    if (v > max) max = v;
    sum += v;
  }
  const mean = sum / data.length;
  let sq = 0;
  for (let i = 0; i < data.length; ++i) sq += (data[i] - mean) ** 2;

  return {
    label,
    grid: {
      transform: { kind: 'matrix', matrix: Mat4.identity() },
      cells: Tensor.create(space, Tensor.Data1(data)),
      stats: { min, max, mean, sigma: Math.sqrt(sq / data.length) },
    },
    instances: [{ transform: Mat4.identity() }],
    sourceData: { kind: 'visual-test', name: label, data: {} } as any,
    customProperties: new CustomProperties(),
    _propertyData: Object.create(null),
    _localPropertyData: Object.create(null),
  };
}

export async function volumeRepr<P extends PD.Params>(
  graphics: GraphicsContext,
  provider: RepresentationProvider<Volume, P>,
  volume: Volume,
  props: Partial<PD.Values<P>>,
  colorTheme?: { name: string; params?: {} },
): Promise<Representation.Any> {
  const registries = getThemeRegistries();
  const repr = provider.factory({ graphics, ...registries }, provider.getParams);
  repr.setTheme(
    Theme.create(
      registries,
      // as the plugin: direct volumes color by value in the shader (`direct-location`), not per cell
      { volume, locationKinds: provider.locationKinds },
      {
        colorTheme: {
          name: colorTheme?.name ?? provider.defaultColorTheme.name,
          params: colorTheme?.params ?? provider.defaultColorTheme.props ?? {},
        },
        sizeTheme: { name: provider.defaultSizeTheme.name, params: provider.defaultSizeTheme.props ?? {} },
      },
    ),
  );
  await repr.createOrUpdate(props, volume).run();
  return repr as unknown as Representation.Any;
}
