/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Text } from '@molstar/graphics/geo/geometry/text/text';
import { Volume } from '@molstar/model/model/volume';
import { DirectVolumeRepresentationProvider } from '@molstar/graphics/repr/volume/direct-volume';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { SliceRepresentationProvider } from '@molstar/graphics/repr/volume/slice';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type CapabilityKey } from '../framework/capabilities.js';
import { type TestView, type VisualTestCase } from '../framework/test-case.js';
import {
  cylinderStar,
  lineCages,
  meshPrimitives,
  meshTorus,
  pointSphere,
  rowTransforms,
  shapeRepr,
  sphereGrid,
  textLabels,
  textureMeshFromMesh,
} from '../scenes/primitives.js';
import { blobVolume, gpuBlobVolume, volumeRepr } from '../scenes/volume.js';
import { helixStructure, structureRepr } from '../scenes/structure.js';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';

export const ObliqueView: TestView = { direction: [-0.4, -0.5, -1] };

/** Sphere and cylinder impostor shaders require frag depth; there is no fallback at the geometry level */
export const Impostor: CapabilityKey[] = ['fragDepth'];

const textProps: Partial<PD.Values<Text.Params>> = {
  fontQuality: 3,
  background: true,
  backgroundOpacity: 0.8,
  backgroundColor: Color(0xdddddd),
  borderWidth: 0.1,
  attachment: 'middle-center',
};

export const GeometryCases: VisualTestCase[] = [
  {
    id: 'geometry/mesh/primitives',
    tags: ['smoke'],
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives())],
  },
  {
    id: 'geometry/mesh/flat-shaded',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives(), { props: { flatShaded: true } })],
  },
  {
    id: 'geometry/mesh/instanced',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshTorus(1), { transforms: rowTransforms(3, 3) })],
  },
  {
    id: 'geometry/spheres/grid',
    tags: ['smoke'],
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(3, 3))],
  },
  {
    id: 'geometry/spheres/instanced',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(2, 2.5), { transforms: rowTransforms(3, 6) })],
  },
  {
    id: 'geometry/spheres/approximate',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(3, 3), { props: { approximate: true } })],
  },
  {
    id: 'geometry/spheres/lod-levels',
    description: 'Per-level `uLod` is set via multi-draw data; everything is in the far level: every 8th sphere at 2x size.',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('spheres', sphereGrid(6, 1.5), {
        props: {
          lodLevels: [
            { minDistance: 0, maxDistance: 1, overlap: 0, stride: 1, scaleBias: 1 },
            { minDistance: 1, maxDistance: 10000000, overlap: 0, stride: 8, scaleBias: 3 },
          ],
        },
      }),
    ],
  },
  {
    id: 'geometry/cylinders/star',
    tags: ['smoke'],
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [await shapeRepr('cylinders', cylinderStar(4), { size: 0.4 })],
  },
  {
    id: 'geometry/points/square',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'square', sizeFactor: 4 } }),
    ],
  },
  {
    id: 'geometry/points/circle',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'circle', sizeFactor: 6 } }),
    ],
  },
  {
    id: 'geometry/points/fuzzy',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'fuzzy', sizeFactor: 8 } }),
    ],
  },
  {
    id: 'geometry/points/attenuation',
    description: 'Point size attenuated by the view distance.',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), {
        props: { pointStyle: 'circle', sizeFactor: 2, pointSizeAttenuation: true },
      }),
    ],
  },
  {
    id: 'geometry/lines/cages',
    tags: ['smoke'],
    view: ObliqueView,
    scene: async () => [await shapeRepr('lines', lineCages(), { props: { sizeFactor: 3 } })],
  },
  {
    id: 'geometry/text/labels',
    tags: ['smoke'],
    description: 'Font atlas is rasterized with a 2D canvas; expect browser/OS dependent snapshots.',
    scene: async () => [await shapeRepr('text', textLabels(textProps), { props: textProps })],
  },
  {
    id: 'geometry/mesh/isosurface-cpu',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, IsosurfaceRepresentationProvider, blobVolume(32), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/primitives',
    description: 'Texture mesh built on the CPU (vertex data in float textures), as GPU-computed surfaces are.',
    view: ObliqueView,
    scene: async ({ graphics }) => [await shapeRepr('texture-mesh', textureMeshFromMesh(graphics, meshPrimitives()))],
  },
  {
    id: 'geometry/texture-mesh/isosurface-gpu',
    tags: ['smoke'],
    requires: ['floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await volumeRepr(webgl, IsosurfaceRepresentationProvider, blobVolume(32), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: true,
      }),
    ],
  },
  {
    id: 'geometry/mesh/grid3d-compute',
    description: 'Volume computed with `createGrid3dComputeRenderable` (loop bound uniform, define on WebGL1), CPU isosurface.',
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await volumeRepr(webgl, IsosurfaceRepresentationProvider, await gpuBlobVolume(webgl, 32, false), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/mesh/grid3d-compute-cumulative',
    description: 'As `grid3d-compute`, accumulated per blob (`tCumulativeSum`).',
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await volumeRepr(webgl, IsosurfaceRepresentationProvider, await gpuBlobVolume(webgl, 32, true), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/mesh/gaussian-surface-cpu',
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(graphics, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
        tryUseGpu: false,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/gaussian-surface-gpu',
    description: 'Gaussian density computed on the GPU (`tMinDistanceTex`, 2d or 3d by `dGridTexType`).',
    requires: ['floatTextureLinear', 'floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await structureRepr(webgl, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
        tryUseGpu: true,
      }),
    ],
  },
  {
    id: 'geometry/texture-mesh/gaussian-surface-gpu-smooth-colors',
    description: 'GPU gaussian surface with GPU color smoothing (`color-accumulate`, `color-normalize`).',
    requires: ['floatTextureLinear', 'floatRenderTarget', 'multipleRenderTargets'],
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await structureRepr(
        webgl,
        GaussianSurfaceRepresentationProvider,
        await helixStructure(24),
        {
          resolution: 0.5,
          tryUseGpu: true,
          smoothColors: { name: 'on', params: { resolutionFactor: 2, sampleStride: 3 } },
        },
        { name: 'sequence-id' },
      ),
    ],
  },
  {
    id: 'geometry/direct-volume/blobs',
    tags: ['smoke'],
    view: ObliqueView,
    scene: async ({ graphics }) => [await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(32), {})],
  },
  {
    id: 'geometry/image/slice',
    scene: async ({ graphics }) => [
      await volumeRepr(graphics, SliceRepresentationProvider, blobVolume(32), {
        dimension: { name: 'relativeZ', params: 0.5 },
      }),
    ],
  },
];
