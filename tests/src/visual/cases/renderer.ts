/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { createClipping } from '@molstar/graphics/geo/geometry/clipping-data';
import type { Representation } from '@molstar/graphics/repr/representation';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { ValueCell } from '@molstar/core/util';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { EveryLoci, Loci } from '@molstar/model/model/loci';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { Mat4, Vec2, Vec3, Vec4 } from '@molstar/core/math/linear-algebra';
import { Clip } from '@molstar/core/util/clip';
import { Color } from '@molstar/core/util/color';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { type SceneContext, type VisualTestCase } from '../framework/test-case.js';
import {
  cylinderStar,
  lineCages,
  meshBackdrop,
  meshBox,
  meshPrimitives,
  meshTorus,
  pointSphere,
  rowTransforms,
  shapeRepr,
  sphereGrid,
} from '../scenes/primitives.js';
import { helixStructure, structureRepr } from '../scenes/structure.js';
import { Impostor, ObliqueView } from './geometry.js';
import { Interval } from '@molstar/core/data/int/interval';
import { OrderedSet } from '@molstar/core/data/int/ordered-set';
import { Volume } from '@molstar/model/model/volume';
import { DirectVolumeRepresentationProvider } from '@molstar/graphics/repr/volume/direct-volume';
import { blobVolume, volumeRepr } from '../scenes/volume.js';

async function primitivesScene() {
  return [await shapeRepr('mesh', meshPrimitives())];
}

async function illuminationScene() {
  return [
    await shapeRepr('mesh', meshPrimitives()),
    await shapeRepr('backdrop', meshBackdrop(12, -2), { colors: [Color(0xcccccc)] }),
  ];
}

/** The scene of `renderer/pick/direct-volume`: a direct volume of blobs and an opaque box beside it */
async function pickScene({ graphics }: SceneContext) {
  const blobs = await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(32), {});
  const box = await shapeRepr('box', meshBox(8), {
    colors: [Color(0x999999)],
    transforms: [Mat4.fromTranslation(Mat4(), Vec3.create(-12, 16, 16))],
  });
  return [blobs, box];
}

/**
 * Picks the pixels of world points (`canvas3d.asyncIdentify`, on every backend): a blob of the volume is a cell where
 * the volume looks solid, the box is the box, the corner nothing; selects the picked volume whole and highlights the box
 */
async function pickDirectVolume({ graphics, canvas3d }: SceneContext, reprs: Representation.Any[]) {
  const [blobs, box] = reprs;
  const errors: string[] = [];
  const pick = async (point: Vec2) => {
    const data = await canvas3d.asyncIdentify(point)?.wait();
    return data ? canvas3d.getLoci(data.id) : undefined;
  };
  const pixelOf = (world: Vec3) => {
    const p = canvas3d.camera.project(Vec4(), world);
    return Vec2.create(p[0] / graphics.pixelRatio, p[1] / graphics.pixelRatio);
  };

  // the center of the third blob of \`blobVolume(32)\`, clear of the box
  const volume = await pick(pixelOf(Vec3.create(16, 16 + 32 / 5, 16 + 32 / 8)));
  if (!volume || !Volume.Cell.isLoci(volume.loci) || volume.repr !== blobs) {
    errors.push(`volume: picked ${volume?.loci.kind ?? 'nothing'}, not a cell of the volume`);
  } else {
    const { volume: v, elements } = volume.loci;
    const cell = OrderedSet.getAt(elements[0].indices, 0);
    const value = v.grid.cells.data[cell];
    if (OrderedSet.size(elements[0].indices) !== 1) errors.push('volume: picked more than one cell');
    if (!(value > v.grid.stats.mean)) errors.push(`volume: picked cell ${cell} of value ${value}, not in a blob`);
    canvas3d.mark({ loci: Loci.applyGranularity(volume.loci, 'structure'), repr: blobs }, MarkerAction.Select);
  }

  const boxPick = await pick(pixelOf(Vec3.create(-12, 16, 16)));
  if (boxPick?.repr !== box) errors.push(`box: picked ${boxPick?.loci.kind ?? 'nothing'}, not the box`);
  else canvas3d.mark(boxPick, MarkerAction.Highlight);

  const corner = await pick(Vec2.create(1, 1));
  if (corner && corner.loci.kind !== 'empty-loci') errors.push(`corner: picked ${corner.loci.kind}, not nothing`);
  return errors;
}

export const RendererCases: VisualTestCase[] = [
  {
    id: 'renderer/context/release-restore',
    description: 'the context released on purpose (a canvas out of view) and restored: the scene is built again on it',
    requires: ['releaseContext', ...Impostor],
    isolate: true,
    releaseContext: true,
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await shapeRepr('mesh', meshPrimitives()),
      await structureRepr(graphics, SpacefillRepresentationProvider, await helixStructure(12), {}),
    ],
  },
  {
    id: 'renderer/context/release-restore-gpu-surface',
    description:
      'the context released and restored with a GPU surface (compute kernels on WebGPU, the vertex arena of the previous device): the surface is computed again on the new context',
    requires: ['releaseContext', 'gaussianSurfaceGpu'],
    isolate: true,
    releaseContext: true,
    view: ObliqueView,
    scene: async ({ graphics }) => [
      await structureRepr(graphics, GaussianSurfaceRepresentationProvider, await helixStructure(24), {
        resolution: 0.5,
        tryUseGpu: true,
      }),
    ],
  },
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
  ...[1, 0.5].map(
    (resolutionScale): VisualTestCase => ({
      id: `renderer/marking/dim-occlusion${resolutionScale === 1 ? '' : '-half-resolution'}`,
      description:
        'Selection fill and dimming of a dense block of spheres behind, shaded with the multi-scale SSAO occlusion' +
        (resolutionScale === 1 ? '.' : ', computed at half resolution and upsampled with the depth.'),
      requires: Impostor,
      view: ObliqueView,
      canvasProps: {
        renderer: { dimStrength: 0.7 },
        postprocessing: {
          occlusion: {
            name: 'on',
            params: {
              ...PD.getDefaultValues(SsaoParams),
              multiScale: { name: 'on', params: SsaoParams.multiScale.map('on').defaultValue },
              resolutionScale,
            },
          },
        },
      },
      scene: async () => {
        const primitives = await shapeRepr('mesh', meshPrimitives());
        primitives.mark(EveryLoci, MarkerAction.Select);
        const behind = Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, -6));
        const block = await shapeRepr('spheres', sphereGrid(6, 1.8), { transforms: [behind] });
        return [primitives, block];
      },
    }),
  ),
  ...[0.3, 0.1].map(
    (surfaceOpacity): VisualTestCase => ({
      id: `renderer/marking/direct-volume${surfaceOpacity === 0.3 ? '' : '-surface-opacity'}`,
      description:
        'A direct volume with its lower half of cells selected around a highlighted opaque box: edges and fill where ' +
        'the volume looks solid (its marking surface), ghost edges behind the box' +
        (surfaceOpacity === 0.3 ? '.' : `; the surface at the opacity ${surfaceOpacity} (\`surfaceOpacity\`): wider.`),
      view: ObliqueView,
      scene: async ({ graphics }) => {
        const n = 32;
        const volume = blobVolume(n);
        const blobs = await volumeRepr(graphics, DirectVolumeRepresentationProvider, volume, { surfaceOpacity });
        const lower = Interval.ofBounds(0, (n * n * n) / 2) as Interval<Volume.CellIndex>;
        blobs.mark(
          Volume.Cell.Loci(volume, [{ indices: lower, instances: Interval.ofSingleton(0 as Volume.InstanceIndex) }]),
          MarkerAction.Select,
        );
        const box = await shapeRepr('box', meshBox(8), {
          colors: [Color(0x999999)],
          transforms: [Mat4.fromTranslation(Mat4(), Vec3.create(16, 16, 16))],
        });
        box.mark(EveryLoci, MarkerAction.Highlight);
        return [blobs, box];
      },
    }),
  ),
  {
    id: 'renderer/pick/direct-volume',
    description:
      'Picking a direct volume where it looks solid (a cell in a blob), the opaque box beside it and the background; ' +
      'the picked volume selected whole, the picked box highlighted.',
    requires: Impostor,
    view: ObliqueView,
    scene: pickScene,
    interact: pickDirectVolume,
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
      'Hi-Z occlusion culling (WebGL2, WebGPU) runs in the interactive frames before the capture; the captured image is not culled.',
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
    id: 'renderer/hi-z/viewport',
    description:
      'Hi-Z occlusion culling in a frame of the canvas (top-left offset); captured from the canvas, so culled objects are missing from the image where the Hi-Z offset is wrong.',
    requires: ['floatRenderTarget'],
    capture: 'canvas',
    canvasProps: {
      hiZ: { enabled: true, minLevel: 1 },
      viewport: { name: 'relative-frame', params: { x: 0.1, y: 0.25, width: 0.75, height: 0.6 } },
    },
    frames: 4,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshTorus(1), { transforms: rowTransforms(5, 3) }),
      await shapeRepr('backdrop', meshBackdrop(12, 2), { colors: [Color(0xcccccc)] }),
    ],
  },
  {
    id: 'renderer/viewport/canvas',
    description: 'Captured from the canvas (Canvas3D frame) instead of an image pass, the whole canvas.',
    capture: 'canvas',
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/viewport/relative-frame',
    description: 'Drawn into a frame relative to the canvas: 10% from the left, 20% from the top (top-left origin).',
    capture: 'canvas',
    canvasProps: { viewport: { name: 'relative-frame', params: { x: 0.1, y: 0.2, width: 0.6, height: 0.5 } } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/viewport/static-frame',
    description: 'Drawn into a frame of 160x112 CSS pixels at 24, 120 from the top-left corner of the canvas.',
    capture: 'canvas',
    canvasProps: { viewport: { name: 'static-frame', params: { x: 24, y: 120, width: 160, height: 112 } } },
    view: ObliqueView,
    scene: primitivesScene,
  },
  {
    id: 'renderer/camera/stereo',
    description: 'Side-by-side stereo: the left eye in the left half, the right eye in the right half.',
    capture: 'canvas',
    size: { width: 384, height: 192 },
    canvasProps: { camera: { stereo: { name: 'on', params: { eyeSeparation: 0.062, focus: 10 } } } },
    view: ObliqueView,
    scene: primitivesScene,
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
    id: 'renderer/clipping/groups',
    description:
      'Per group clipping flags exempting groups from clip objects, from the third object on: the sphere is exempt ' +
      'from the third object (a plane), the box from the first and second, the torus from the second and third.',
    view: ObliqueView,
    scene: async () => {
      const repr = await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel', true) } });
      const { One, Two, Three, None } = Clipping.Groups.Flag;
      setClippingGroups(repr, [Three, One | Two, Two | Three, None]);
      return [repr];
    },
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

/** Clipping flags per group of the single instance of `repr`, as `Visual.setClipping` writes them */
function setClippingGroups(repr: Representation.Any, flags: Clipping.Groups[]) {
  const { values } = repr.renderObjects[0];
  createClipping(flags.length, 'groupInstance', values);
  values.tClipping.ref.value.array.set(flags);
  ValueCell.update(values.tClipping, values.tClipping.ref.value);
}

/** a tilted plane and a sphere (and a horizontal plane with `third`), exercises clip object arrays */
function clipProps(variant: Clip.Variant, third = false): Clip.Props {
  const objects: Clip.Props['objects'] = [
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
  ];
  if (third) {
    objects.push({
      type: 'plane',
      invert: false,
      position: Vec3.create(0, 1, 0),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 90 },
      scale: Vec3.create(1, 1, 1),
      transform: Mat4.identity(),
    });
  }
  return { variant, objects };
}
