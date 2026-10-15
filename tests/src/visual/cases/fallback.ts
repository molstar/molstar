/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { IsosurfaceRepresentationProvider } from '@molstar/graphics/repr/volume/isosurface';
import { Color } from '@molstar/core/util/color';
import { Volume } from '@molstar/model/model/volume';
import { type VisualTestCase } from '../framework/test-case.js';
import { meshBackdrop, meshBox, meshPrimitives, meshTorus, shapeRepr } from '../scenes/primitives.js';
import { helixStructure, structureRepr } from '../scenes/structure.js';
import { blobVolume, volumeRepr } from '../scenes/volume.js';
import { ObliqueView } from './geometry.js';

/**
 * Features requested on profiles without the needed capabilities must fall back to the
 * documented alternative
 */

async function transparentMeshes() {
  return [
    await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] }),
    await shapeRepr('torus', meshTorus(3), { colors: [Color(0x377eb8)], props: { alpha: 0.6 } }),
  ];
}

export const FallbackCases: VisualTestCase[] = [
  {
    id: 'fallback/transparency/wboit',
    contextProps: { transparency: 'wboit' },
    fallback: { requires: ['wboit'], overrides: { contextProps: { transparency: 'blended' } } },
    view: ObliqueView,
    scene: transparentMeshes,
  },
  {
    id: 'fallback/transparency/dpoit',
    contextProps: { transparency: 'dpoit' },
    fallback: { requires: ['dpoit'], overrides: { contextProps: { transparency: 'blended' } } },
    view: ObliqueView,
    scene: transparentMeshes,
  },
  {
    id: 'fallback/illumination',
    canvasProps: { illumination: { enabled: true, maxIterations: 3 } },
    fallback: { requires: ['illumination'], overrides: { canvasProps: { illumination: { enabled: false } } } },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives()),
      await shapeRepr('backdrop', meshBackdrop(12, -2), { colors: [Color(0xcccccc)] }),
    ],
  },
  {
    id: 'fallback/isosurface-gpu',
    fallback: {
      requires: ['floatRenderTarget', 'multipleRenderTargets'],
      overrides: {
        scene: async ({ webgl }) => [
          await volumeRepr(webgl, IsosurfaceRepresentationProvider, blobVolume(32), {
            isoValue: Volume.IsoValue.absolute(0.5),
            tryUseGpu: false,
          }),
        ],
      },
    },
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await volumeRepr(webgl, IsosurfaceRepresentationProvider, blobVolume(32), {
        isoValue: Volume.IsoValue.absolute(0.5),
        tryUseGpu: true,
      }),
    ],
  },
  {
    id: 'fallback/gaussian-surface-gpu',
    fallback: {
      requires: ['floatTextureLinear', 'floatRenderTarget', 'multipleRenderTargets'],
      overrides: {
        scene: async ({ webgl }) => [
          await structureRepr(webgl, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
            resolution: 0.5,
            tryUseGpu: false,
          }),
        ],
      },
    },
    view: ObliqueView,
    scene: async ({ webgl }) => [
      await structureRepr(webgl, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
        tryUseGpu: true,
      }),
    ],
  },
];
