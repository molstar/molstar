/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { EveryLoci } from '@molstar/model/model/loci';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Clip } from '@molstar/core/util/clip';
import { Color } from '@molstar/core/util/color';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { type VisualTestCase } from '../framework/test-case.js';
import {
  cylinderStar,
  lineCages,
  meshBackdrop,
  meshPrimitives,
  meshTorus,
  pointSphere,
  rowTransforms,
  shapeRepr,
  sphereGrid,
} from '../scenes/primitives.js';
import { Impostor, ObliqueView } from './geometry.js';

async function primitivesScene() {
  return [await shapeRepr('mesh', meshPrimitives())];
}

async function illuminationScene() {
  return [
    await shapeRepr('mesh', meshPrimitives()),
    await shapeRepr('backdrop', meshBackdrop(12, -2), { colors: [Color(0xcccccc)] }),
  ];
}

export const RendererCases: VisualTestCase[] = [
  {
    id: 'renderer/camera/orthographic',
    canvasProps: { camera: { mode: 'orthographic' } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/camera/fov-90',
    canvasProps: { camera: { fov: 90 } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/camera/fog',
    tags: ['smoke'],
    requires: Impostor,
    canvasProps: { cameraFog: { name: 'on', params: { intensity: 50 } } },
    view: { direction: [0, -1, -0.3] },
    scene: async () => [await shapeRepr('spheres', sphereGrid(4, 3))],
  },
  {
    id: 'renderer/background/black',
    canvasProps: { renderer: { backgroundColor: Color(0x000000) } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/lighting/two-lights',
    canvasProps: {
      renderer: {
        light: [
          { inclination: 150, azimuth: 320, color: Color(0xffdddd), intensity: 0.5 },
          { inclination: 60, azimuth: 120, color: Color(0xddddff), intensity: 0.4 },
        ],
      },
    },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/lighting/exposure',
    canvasProps: { renderer: { exposure: 1.6, ambientIntensity: 0.6 } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/shading/ignore-light',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives(), { props: { ignoreLight: true } })],
  },
  {
    id: 'renderer/shading/cel',
    view: ObliqueView,
    canvasProps: { renderer: { celSteps: 3 } },
    scene: async () => [await shapeRepr('mesh', meshPrimitives(), { props: { celShaded: true } })],
  },
  {
    id: 'renderer/shading/xray',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives(), { props: { xrayShaded: true } })],
  },
  {
    id: 'renderer/shading/xray-inverted',
    view: ObliqueView,
    scene: async () => [await shapeRepr('mesh', meshPrimitives(), { props: { xrayShaded: 'inverted' } })],
  },
  {
    id: 'renderer/marking/highlight-select',
    tags: ['smoke'],
    requires: Impostor,
    view: ObliqueView,
    scene: async () => {
      const highlighted = await shapeRepr('mesh', meshPrimitives());
      highlighted.mark(EveryLoci, MarkerAction.Highlight);
      const selected = await shapeRepr('spheres', sphereGrid(2, 2), { props: { alpha: 1 } });
      selected.mark(EveryLoci, MarkerAction.Select);
      return [highlighted, selected];
    },
  },
  {
    id: 'renderer/marking/alpha-thickness',
    description:
      'Selected transparent spheres whose opacity scales with their thickness (`alphaThickness`), as does the ' +
      'opacity of their selection fill.',
    requires: Impostor,
    // the fill of the faint spheres differs by a few levels when the opacity is wrong
    tolerance: { threshold: 0.02 },
    view: ObliqueView,
    scene: async () => {
      const selected = await shapeRepr('spheres', sphereGrid(3, 2.5), { props: { alpha: 0.6, alphaThickness: 4 } });
      selected.mark(EveryLoci, MarkerAction.Select);
      return [selected];
    },
  },
  {
    id: 'renderer/multisample/on',
    canvasProps: { multiSample: { mode: 'on', sampleLevel: 2 } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/illumination/default',
    tags: ['smoke'],
    requires: ['illumination'],
    canvasProps: { illumination: { enabled: true, maxIterations: 3 } },
    view: ObliqueView,
    scene: illuminationScene,
  },
  {
    id: 'renderer/illumination/shadows',
    requires: ['illumination'],
    canvasProps: { illumination: { enabled: true, maxIterations: 3, shadowEnable: true } },
    view: ObliqueView,
    scene: illuminationScene,
  },
  {
    id: 'renderer/illumination/primitives',
    description: 'Spheres, points and lines in the traced image.',
    requires: ['illumination', ...Impostor],
    canvasProps: { illumination: { enabled: true, maxIterations: 3 } },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('spheres', sphereGrid(2, 3)),
      await shapeRepr('points', pointSphere(400, 4), { props: { pointStyle: 'circle', sizeFactor: 4 } }),
      await shapeRepr('lines', lineCages()),
      await shapeRepr('backdrop', meshBackdrop(12, -4), { colors: [Color(0xcccccc)] }),
    ],
  },
  {
    id: 'renderer/illumination/transparency',
    description: 'Transparent spheres composed over the traced opaque scene.',
    requires: ['illumination', ...Impostor],
    canvasProps: { illumination: { enabled: true, maxIterations: 3 } },
    view: ObliqueView,
    scene: async () => [
      ...(await illuminationScene()),
      await shapeRepr('spheres', sphereGrid(2, 3), { props: { alpha: 0.5 } }),
    ],
  },
  {
    id: 'renderer/illumination/marking',
    description: 'Selection fill and dimming with the traced shading.',
    requires: ['illumination'],
    canvasProps: { illumination: { enabled: true, maxIterations: 3 }, renderer: { dimStrength: 0.5 } },
    view: ObliqueView,
    scene: async () => {
      const [primitives, backdrop] = await illuminationScene();
      primitives.mark(EveryLoci, MarkerAction.Select);
      return [primitives, backdrop];
    },
  },
  {
    id: 'renderer/illumination/multisample',
    requires: ['illumination'],
    canvasProps: { illumination: { enabled: true, maxIterations: 3 }, multiSample: { mode: 'on', sampleLevel: 2 } },
    view: ObliqueView,
    scene: illuminationScene,
  },
  {
    id: 'renderer/illumination/orthographic-transparent-background',
    requires: ['illumination'],
    canvasProps: {
      illumination: { enabled: true, maxIterations: 3 },
      camera: { mode: 'orthographic' },
      transparentBackground: true,
    },
    view: ObliqueView,
    scene: illuminationScene,
  },
  {
    id: 'renderer/illumination/occlusion-outline',
    description: 'Occlusion and outlines (not ignored) composed with the traced image and transparent spheres.',
    requires: ['illumination', ...Impostor],
    canvasProps: {
      illumination: { enabled: true, maxIterations: 3, ignoreOutline: false },
      postprocessing: {
        occlusion: { name: 'on', params: PD.getDefaultValues(SsaoParams) },
        outline: { name: 'on', params: PD.getDefaultValues(OutlineParams) },
      },
    },
    view: ObliqueView,
    scene: async () => [
      ...(await illuminationScene()),
      await shapeRepr('spheres', sphereGrid(2, 3), { props: { alpha: 0.5 } }),
    ],
  },
  {
    id: 'renderer/hi-z/occlusion',
    description:
      'Hi-Z occlusion culling (WebGL2 only) runs in the interactive frames before the capture; the captured image is not culled.',
    requires: ['floatRenderTarget'],
    canvasProps: { hiZ: { enabled: true, minLevel: 1 } },
    frames: 4,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshTorus(1), { transforms: rowTransforms(5, 3) }),
      await shapeRepr('backdrop', meshBackdrop(12, 2), { colors: [Color(0xcccccc)] }),
    ],
  },
  {
    id: 'renderer/clipping/pixel',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel') } }),
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clipping/solid-interior',
    description: 'Opaque caps of clipped meshes, stencil-marked and filled by the renderer.',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel'), solidInterior: true } }),
    ],
  },
  {
    id: 'renderer/clipping/impostor-interior',
    description: 'Solid interiors of clipped sphere and cylinder impostors, capped on the clip object surfaces.',
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel'), solidInterior: true } }),
      await shapeRepr('cylinders', cylinderStar(4), {
        size: 0.4,
        props: { clip: clipProps('pixel'), solidInterior: true },
      }),
    ],
  },
  {
    id: 'renderer/clipping/points-lines',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), {
        props: { pointStyle: 'circle', sizeFactor: 6, clip: clipProps('pixel') },
      }),
      await shapeRepr('lines', lineCages(), { props: { sizeFactor: 3, clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clipping/instance',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshTorus(1), {
        transforms: rowTransforms(5, 3),
        props: { clip: clipProps('instance') },
      }),
    ],
  },
];

/** a tilted plane and a sphere, exercises clip object arrays */
function clipProps(variant: Clip.Variant): Clip.Props {
  return {
    variant,
    objects: [
      {
        type: 'plane',
        invert: false,
        position: Vec3.create(0, 0, 0),
        rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
        scale: Vec3.create(1, 1, 1),
        transform: Mat4.identity(),
      },
      {
        type: 'sphere',
        invert: false,
        position: Vec3.create(-2, -2, 0),
        rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
        scale: Vec3.create(3, 3, 3),
        transform: Mat4.identity(),
      },
    ],
  };
}
