/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Mat4, Tensor, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import { Volume } from '@molstar/model/model/volume';
import { DirectVolumeRepresentationProvider } from '@molstar/graphics/repr/volume/direct-volume';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { GaussianVolumeRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-volume';
import { blobVolume, createVolume, volumeRepr } from '../visual/scenes/volume.js';
import { structureRepr } from '../visual/scenes/structure.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { MarkingParams } from '@molstar/graphics/canvas3d/passes/marking';
import type { BenchmarkCase } from './framework.js';
import { LargeSize, type RenderInput, renderCase } from './render.js';
import { exampleStructure } from './structures.js';

/** deterministic uniform random numbers in [0, 1) (mulberry32) */
function random(seed: number) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

/**
 * A map-like volume on a `dim` grid: `clusters` globules of `atoms` Gaussian atoms each (sigma 1.2 cells, globule
 * radius 12 cells) over uniform noise; most cells are background, as in cryo-EM maps and tomograms. `matrix` is the
 * grid's transform (cell axes and spacing)
 */
function mapVolume(label: string, dim: Vec3, clusters: number, atoms: number, matrix = Mat4.identity()): Volume {
  const [nx, ny, nz] = dim;
  const space = Tensor.Space([nx, ny, nz], [2, 1, 0], Float32Array);
  const data = new Float32Array(nx * ny * nz);
  const rand = random(42);
  for (let i = 0; i < data.length; ++i) data[i] = rand() * 0.05;

  const sigma = 1.2;
  const r = Math.ceil(3 * sigma);
  const f = 1 / (2 * sigma * sigma);
  const radius = 12;
  const margin = radius + r;
  for (let c = 0; c < clusters; ++c) {
    const cx = margin + rand() * (nx - 2 * margin);
    const cy = margin + rand() * (ny - 2 * margin);
    const cz = Math.min(nz / 2, margin) + rand() * Math.max(0, nz - 2 * margin);
    for (let a = 0; a < atoms; ++a) {
      // uniform in a ball
      let dx: number, dy: number, dz: number;
      do {
        dx = rand() * 2 - 1;
        dy = rand() * 2 - 1;
        dz = rand() * 2 - 1;
      } while (dx * dx + dy * dy + dz * dz > 1);
      const x = cx + dx * radius;
      const y = cy + dy * radius;
      const z = cz + dz * radius;
      const x0 = Math.max(0, Math.floor(x - r)),
        x1 = Math.min(nx - 1, Math.ceil(x + r));
      const y0 = Math.max(0, Math.floor(y - r)),
        y1 = Math.min(ny - 1, Math.ceil(y + r));
      const z0 = Math.max(0, Math.floor(z - r)),
        z1 = Math.min(nz - 1, Math.ceil(z + r));
      for (let i = x0; i <= x1; ++i) {
        for (let j = y0; j <= y1; ++j) {
          for (let k = z0; k <= z1; ++k) {
            const d2 = (i - x) ** 2 + (j - y) ** 2 + (k - z) ** 2;
            data[space.dataOffset(i, j, k)] += Math.exp(-d2 * f);
          }
        }
      }
    }
  }
  const volume = createVolume(label, space, data);
  return { ...volume, grid: { ...volume.grid, transform: { kind: 'matrix', matrix } } };
}

const cached = <T>(create: () => T) => {
  let value: T | undefined;
  return () => (value ??= create());
};

const getBlobs = cached(() => blobVolume(256));
const getMap = cached(() => mapVolume('map-256', Vec3.create(256, 256, 256), 40, 300));
/** a slab, tilted by 30° about x: the grid's bounding sphere far larger than its depth */
const getTomogram = cached(() =>
  mapVolume(
    'tomogram-512x512x96',
    Vec3.create(512, 512, 96),
    150,
    300,
    Mat4.fromRotation(Mat4(), Math.PI / 6, Vec3.create(1, 0, 0)),
  ),
);
/** the cube of `getMap` sampled with cells of 1x1x4 */
const getAnisotropic = cached(() =>
  mapVolume('map-256x256x64', Vec3.create(256, 256, 64), 40, 300, Mat4.fromScaling(Mat4(), Vec3.create(1, 1, 4))),
);

/** a wide, faint transfer function: most cells of the blobs contribute and rays rarely saturate */
const DenseControlPoints = [Vec2.create(0.05, 0), Vec2.create(0.1, 0.01), Vec2.create(0.99, 0.01)];

function volumeCase(
  id: string,
  description: string,
  volume: () => Volume,
  props: object = {},
  input: Partial<RenderInput> = {},
) {
  return renderCase(`render/direct-volume/${id}`, `${description}, 2560x1440`, async () => {
    const v = volume();
    return {
      size: LargeSize,
      scene: async (g) => [await volumeRepr(g, DirectVolumeRepresentationProvider, v, props)],
      ...input,
    };
  });
}

/**
 * Direct-volume rendering at 2560x1440 (`--gpu-timing`: the `volume` pass), not in the default runs
 * (`--filter 'render/direct-volume/*'`)
 */
export const DirectVolumeCases: BenchmarkCase<RenderInput>[] = [
  volumeCase('blobs-256-dense', 'Three Gaussian blobs on a 256³ grid, a wide faint transfer function', getBlobs, {
    controlPoints: DenseControlPoints,
  }),
  volumeCase('map-256', 'A map-like 256³ grid (40 globules over noise), the default ramp', getMap),
  volumeCase(
    'map-256-zoom',
    'As map-256, the camera at 0.4 of its distance (the volume beyond the screen)',
    getMap,
    {},
    {
      zoom: 0.4,
    },
  ),
  volumeCase('map-256-float', 'As map-256, float data (16 bytes per cell)', getMap, { dataType: 'float' }),
  volumeCase(
    'map-256-highlight',
    'As map-256, highlighted with the marking pass off (markers colored in the volume shader)',
    getMap,
    {},
    { highlight: true, props: { marking: { ...PD.getDefaultValues(MarkingParams), enabled: false } } },
  ),
  volumeCase(
    'map-256-highlight-marking',
    'As map-256, highlighted with the marking pass on (the mask drawn from the marking surface of the volume)',
    getMap,
    {},
    { highlight: true },
  ),
  renderCase(
    'render/direct-volume/map-256-isosurface',
    'As map-256 with an opaque isosurface of the map inside (the volume at the edges of opaque objects), 2560x1440',
    async () => {
      const v = getMap();
      const { mean, sigma } = v.grid.stats;
      return {
        size: LargeSize,
        scene: async (g) => [
          await volumeRepr(g, IsosurfaceRepresentationProvider, v, {
            isoValue: Volume.IsoValue.absolute(mean + 6 * sigma),
          }),
          await volumeRepr(g, DirectVolumeRepresentationProvider, v, {}),
        ],
      };
    },
  ),
  volumeCase('tomogram-512x512x96', 'A map-like slab of 512x512x96 cells (150 globules), tilted by 30°', getTomogram),
  volumeCase(
    'map-256x256x64-anisotropic',
    'A map-like grid of 256x256x64 cells of 1x1x4 (40 globules), the default ramp',
    getAnisotropic,
  ),
  renderCase(
    'render/direct-volume/7qpd-gaussian-volume',
    'PDB entry 7qpd as Gaussian volume (GPU density, packed groups, gradients per sample), 2560x1440',
    async () => {
      const s = await exampleStructure('7qpd.fw2.cif');
      return {
        size: LargeSize,
        scene: async (g) => [await structureRepr(g, GaussianVolumeRepresentationProvider, s, {})],
      };
    },
  ),
];
