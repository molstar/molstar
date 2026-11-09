/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3D, Canvas3DParams } from '@molstar/graphics/canvas3d/canvas3d';
import { now } from '@molstar/core/util/now';
import { DefaultAssetManager } from '@molstar/core/util/assets';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { createClipping } from '@molstar/graphics/geo/geometry/clipping-data';
import type { Representation } from '@molstar/graphics/repr/representation';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { ValueCell } from '@molstar/core/util';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { MarkingParams } from '@molstar/graphics/canvas3d/passes/marking';
import { ClipGhostParams } from '@molstar/graphics/canvas3d/passes/clip-ghost';
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
  cylinderGrid,
  cylinderStar,
  lineCages,
  lineStrips,
  meshBackdrop,
  meshBox,
  meshPrimitives,
  meshTorus,
  gridTransforms,
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

/** The ghost of what clip objects remove, its look by default */
const ClipGhostOn = { name: 'on' as const, params: PD.getDefaultValues(ClipGhostParams) };

export const RendererCases: VisualTestCase[] = [
  {
    id: 'renderer/context/shared-canvas',
    description:
      'Two Canvas3D instances of one context, each with half of its canvas: a frame of one keeps the image of the other.',
    scene: async () => [],
    check: async (_, profile) => {
      const canvas = document.createElement('canvas');
      canvas.width = canvas.height = 256;
      const container = document.createElement('div');
      container.style.cssText = 'position: absolute; left: -10000px; width: 256px; height: 256px';
      container.appendChild(canvas);
      document.body.appendChild(container);
      const context = await profile.createContext(canvas, new DefaultAssetManager());
      if (!context) return [`no context of profile '${profile.id}'`];
      const create = async (x: number, backgroundColor: Color) => {
        const c = Canvas3D.create(context, {
          viewport: { name: 'relative-frame', params: { x, y: 0, width: 0.5, height: 1 } },
          renderer: { ...PD.getDefaultValues(Canvas3DParams).renderer, backgroundColor },
        });
        c.add(await shapeRepr('mesh', meshPrimitives()));
        c.commit(true);
        return c;
      };
      const frame = async (c: Canvas3D) => {
        await context.graphics.resources.pending.whenReady();
        c.requestDraw();
        c.tick(now(), { isSynchronous: true });
      };
      const read = () => {
        const copy = document.createElement('canvas');
        copy.width = copy.height = 256;
        const ctx = copy.getContext('2d')!;
        ctx.drawImage(canvas, 0, 0);
        const at = (x: number) => Array.from(ctx.getImageData(x, 2, 1, 1).data).join(',');
        return { left: at(2), right: at(253) };
      };
      const left = await create(0, Color(0xff0000));
      // what an instance drew alone is not kept: one added asks the others for a frame
      let requested = 0;
      const removeProbe = context.addInstance!(() => {
        requested += 1;
      });
      const right = await create(0.5, Color(0x0000ff));
      removeProbe();
      try {
        for (let i = 0; i < 3; ++i) {
          await frame(left);
          await frame(right);
          await context.graphics.waitForGpuCommandsComplete();
          await new Promise((r) => requestAnimationFrame(r));
        }
        await frame(left);
        await frame(right);
        const both = read();
        await context.graphics.waitForGpuCommandsComplete();
        await new Promise((r) => requestAnimationFrame(r));
        // a new frame, of one instance only
        await frame(left);
        const one = read();
        const red = '255,0,0,255';
        const blue = '0,0,255,255';
        const errors: string[] = [];
        if (requested !== 1) errors.push(`frames asked of an instance when another is added: ${requested}`);
        if (both.left !== red || both.right !== blue) errors.push(`both drawn: ${both.left} | ${both.right}`);
        if (one.left !== red || one.right !== blue) errors.push(`one drawn: ${one.left} | ${one.right}`);
        // a resize leaves nothing to keep; frames go on without the other instance
        canvas.width = canvas.height = 128;
        await context.graphics.waitForGpuCommandsComplete();
        await new Promise((r) => requestAnimationFrame(r));
        await frame(left);
        await context.graphics.waitForGpuCommandsComplete();
        right.dispose();
        await new Promise((r) => requestAnimationFrame(r));
        await frame(left);
        await context.graphics.waitForGpuCommandsComplete();
        errors.push(...((context.graphics as { takeErrors?: () => string[] }).takeErrors?.() ?? []));
        return errors;
      } finally {
        left.dispose();
        right.dispose();
        context.dispose?.();
        container.remove();
      }
    },
  },
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
  ...[true, false].map(
    (reuseOcclusion): VisualTestCase => ({
      id: `renderer/multisample/ssao${reuseOcclusion ? '-reused' : ''}`,
      description:
        'A block of spheres with SSAO, multi-sampled (8 samples): ' +
        (reuseOcclusion
          ? "the first sample's occlusion read by the others where each one's pixels lie (`reuseOcclusion`), " +
            'closer to the occlusion of each sample than without the shift or with its y the other way.'
          : 'the occlusion of each sample.'),
      requires: Impostor,
      view: ObliqueView,
      canvasProps: {
        multiSample: { mode: 'on', sampleLevel: 3, reuseOcclusion },
        postprocessing: {
          occlusion: { name: 'on', params: { ...PD.getDefaultValues(SsaoParams), resolutionScale: 1 } },
        },
      },
      scene: async () => [await shapeRepr('spheres', sphereGrid(6, 1.8))],
    }),
  ),
  ...(['on', 'temporal'] as const).map(
    (mode): VisualTestCase => ({
      id: `renderer/multisample/canvas-${mode}`,
      description:
        mode === 'on'
          ? 'Mesh primitives and spheres multi-sampled (4 samples), the canvas captured: the samples composed into the drawing buffer.'
          : 'Mesh primitives and spheres with temporal multi-sampling, the canvas captured after a forced frame (as while ' +
            "the camera moves): the unjittered first sample presented from the draw pass's color target.",
      requires: Impostor,
      capture: 'canvas',
      view: ObliqueView,
      canvasProps: { multiSample: { mode, sampleLevel: 2 } },
      scene: primitivesScene,
    }),
  ),
  {
    id: 'renderer/marking/multi-sample',
    description:
      'Highlighted mesh primitives and selected spheres, multi-sampled: the marking of each jittered sample tested ' +
      "against the scene's depth drawn with another jitter, no ghost slices on the marked surfaces.",
    requires: Impostor,
    view: ObliqueView,
    canvasProps: { multiSample: { mode: 'on', sampleLevel: 2 } },
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
  ...[1, 0.5].flatMap((resolutionScale): VisualTestCase[] => {
    const suffix = resolutionScale === 1 ? '' : '-half-resolution';
    const canvasProps: VisualTestCase['canvasProps'] = {
      renderer: { dimStrength: 0.7 },
      postprocessing: {
        occlusion: {
          name: 'on',
          params: { ...PD.getDefaultValues(SsaoParams), transparentThreshold: 1, resolutionScale },
        },
      },
      // the layer (alpha 0.5) hides the marked block behind it: read from the transparent depth
      marking: { ...PD.getDefaultValues(MarkingParams), transparentThreshold: 1 },
    };
    // a transparent layer of spheres in front of an opaque block
    const layerScene = async (mark: 'layer' | 'block') => {
      const front = Mat4.fromTranslation(Mat4(), Vec3.create(0, 0, 6));
      const layer = await shapeRepr('spheres', sphereGrid(3, 2.5), { transforms: [front], props: { alpha: 0.5 } });
      const block = await shapeRepr('spheres', sphereGrid(6, 1.8));
      (mark === 'layer' ? layer : block).mark(EveryLoci, MarkerAction.Select);
      return [layer, block];
    };
    return [
      {
        id: `renderer/marking/transparent-layer-marked${suffix}`,
        description:
          'A selected transparent layer over a dimmed opaque block: its fill takes the occlusion of the layer, not ' +
          'the darker one of the block seen through it.',
        requires: Impostor,
        view: ObliqueView,
        canvasProps,
        scene: () => layerScene('layer'),
      },
      {
        id: `renderer/marking/transparent-layer-unmarked${suffix}`,
        description:
          'A selected opaque block behind an unmarked transparent layer: its fill is attenuated by the layer, which ' +
          'is itself dimmed by its alpha.',
        requires: Impostor,
        view: ObliqueView,
        canvasProps,
        scene: () => layerScene('block'),
      },
    ];
  }),
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
    id: 'renderer/viewport/outside-untouched',
    description:
      'A frame relative to the canvas: nothing is drawn outside it, which a context may share with other Canvas3D instances.',
    // a canvas of its own: another case's frame stays outside the frame on WebGL (`preserveDrawingBuffer`)
    isolate: true,
    capture: 'canvas',
    canvasProps: {
      renderer: { backgroundColor: Color(0xff00ff) },
      viewport: { name: 'relative-frame', params: { x: 0.1, y: 0.2, width: 0.6, height: 0.6 } },
    },
    checkImage: ({ data }) => {
      const corner = Array.from(data.subarray(0, 4)).join(',');
      return corner === '0,0,0,0' ? [] : [`drawn outside the frame: ${corner}`];
    },
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
    id: 'renderer/camera/stereo-spheres',
    description: 'Side-by-side stereo with a wide eye separation: sphere impostors in strongly off-axis projections.',
    requires: Impostor,
    capture: 'canvas',
    size: { width: 384, height: 192 },
    canvasProps: { camera: { stereo: { name: 'on', params: { eyeSeparation: 1, focus: 3 } } } },
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(3, 3))],
  },
  ...(
    [
      ['default', {}, [1.25, 0.6, 3]],
      ['solid-interior', { solidInterior: true }, [1.25, 0.6, 3]],
      ['double-sided', { solidInterior: false, doubleSided: true }, [1.25, 0.6, 3]],
      ['open', { solidInterior: false }, [1.25, 0.6, 3]],
      ['eye-inside', {}, [2.6, 0.2, 2.7]],
      ['eye-inside-open', { solidInterior: false }, [2.6, 0.2, 2.7]],
    ] as const
  ).map(([name, props, [x, y, z]]) => ({
    id: `renderer/sphere-near/${name}`,
    description:
      'The camera in a grid of sphere impostors: spheres cut by the near plane and reaching behind the eye, their ' +
      `quads the bounds of the part beyond the plane (${name}).`,
    requires: Impostor,
    camera: {
      position: Vec3.create(x, y, z),
      target: Vec3.create(x, y, -5),
      up: Vec3.create(0, 1, 0),
      radius: 3,
      radiusMax: 3,
    },
    scene: async () => [await shapeRepr('spheres', sphereGrid(5, 2.5), { size: 1.2, props })],
  })),
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
  ...[false, true].map((solidInterior) => ({
    id: `renderer/clipping/impostor-classes${solidInterior ? '-interior' : ''}`,
    description:
      'Sphere and cylinder impostors on a grid by pixel clipping: clipped whole (moved out in the vertex stage), cut ' +
      `(tested per fragment) and not cut (not tested)${solidInterior ? ', with solid interiors' : ''}.`,
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('spheres', sphereGrid(4, 2), { size: 0.5, props: { clip: clipProps('pixel'), solidInterior } }),
      await shapeRepr('cylinders', cylinderGrid(3, 3, 1.6), {
        size: 0.3,
        props: { clip: clipProps('pixel'), solidInterior },
      }),
    ],
  })),
  {
    id: 'renderer/clipping/points-lines',
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), {
        props: { pointStyle: 'circle', sizeFactor: 6, clip: clipProps('pixel') },
      }),
      await shapeRepr('lines', lineCages(), { props: { sizeFactor: 3, clip: clipProps('pixel') } }),
      await shapeRepr('lines', lineStrips(), { props: { sizeFactor: 3, clip: clipProps('pixel') } }),
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
  ...(['plane', 'inverted-sphere'] as const).map((clip) => ({
    id: `renderer/clipping/solid-interior-cells-${clip}`,
    description:
      'Solid interiors of instanced tori in small instance-grid cells, those clipped whole culled: the caps of the cut ' +
      `ones counted over every face behind them${clip === 'inverted-sphere' ? ' (the removed outside of a kept sphere, behind it too)' : ''}.`,
    requires: Impostor,
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshTorus(0.8), {
        transforms: gridTransforms(5, 2.2),
        props: {
          solidInterior: true,
          cellSize: 2,
          batchSize: 6,
          clip: {
            variant: 'pixel' as const,
            objects: [
              clip === 'plane'
                ? {
                    type: 'plane' as const,
                    invert: false,
                    position: Vec3.create(0.3, 0, 0),
                    rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
                    scale: Vec3.create(1, 1, 1),
                    transform: Mat4.identity(),
                  }
                : {
                    type: 'sphere' as const,
                    invert: true,
                    position: Vec3.create(0.4, 0.2, 0),
                    rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
                    scale: Vec3.create(7, 7, 7),
                    transform: Mat4.identity(),
                  },
            ],
          },
        },
      }),
    ],
  })),
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
  {
    id: 'renderer/clip-ghost/pixel',
    description:
      'What the clip objects remove as a ghost (a tint and its edges at half resolution): of a mesh by fragment, of ' +
      'sphere impostors whole by their centers.',
    requires: Impostor,
    canvasProps: { clipGhost: ClipGhostOn },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel') } }),
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clip-ghost/points-lines',
    description: 'The ghost of points and lines, each whole by its center.',
    canvasProps: { clipGhost: ClipGhostOn },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('points', pointSphere(400, 4), {
        props: { pointStyle: 'circle', sizeFactor: 6, clip: clipProps('pixel') },
      }),
      await shapeRepr('lines', lineCages(), { props: { sizeFactor: 3, clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clip-ghost/instance',
    description: 'The ghost of instances clipped by their centers: the removed tori whole.',
    canvasProps: { clipGhost: ClipGhostOn },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshTorus(1), {
        transforms: rowTransforms(5, 3),
        props: { clip: clipProps('instance') },
      }),
    ],
  },
  {
    id: 'renderer/clip-ghost/grid-cells',
    description:
      'The ghost of instanced tori in small instance-grid cells cut by a plane: the cells clipped whole culled from ' +
      'the visible pass, drawn by the ghost pass (and the cells the plane does not reach the other way round).',
    canvasProps: { clipGhost: ClipGhostOn },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshTorus(0.8), {
        transforms: gridTransforms(5, 2.2),
        props: {
          cellSize: 2,
          batchSize: 6,
          clip: {
            variant: 'pixel' as const,
            objects: [
              {
                type: 'plane' as const,
                invert: false,
                position: Vec3.create(0.3, 0, 0),
                rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
                scale: Vec3.create(1, 1, 1),
                transform: Mat4.identity(),
              },
            ],
          },
        },
      }),
    ],
  },
  {
    id: 'renderer/clip-ghost/opt-out',
    description: 'A mesh opted out of the ghost (`clipGhost` false): only the spheres have one.',
    requires: Impostor,
    canvasProps: { clipGhost: ClipGhostOn },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel'), clipGhost: false } }),
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clip-ghost/illumination',
    description: 'The ghost over the traced image, as over the visible passes’.',
    requires: ['illumination', ...Impostor],
    canvasProps: { clipGhost: ClipGhostOn, illumination: { enabled: true, maxIterations: 3 } },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel') } }),
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clip-ghost/illumination-transparent',
    description:
      'The ghost over the traced image with a transparent mesh: behind it where it is fainter than the marking ' +
      'threshold, hidden where denser (the transparent depth written for it).',
    requires: ['illumination', ...Impostor],
    canvasProps: { clipGhost: ClipGhostOn, illumination: { enabled: true, maxIterations: 3 } },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel'), alpha: 0.3 } }),
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel') } }),
    ],
  },
  {
    id: 'renderer/clip-ghost/full-resolution',
    description: 'The ghost at the full resolution, darker: its edges a pixel wide.',
    requires: Impostor,
    canvasProps: {
      clipGhost: { name: 'on', params: { color: Color(0x404040), alpha: 0.3, edgeAlpha: 0.8, resolution: 1 } },
    },
    view: ObliqueView,
    scene: async () => [
      await shapeRepr('mesh', meshPrimitives(), { props: { clip: clipProps('pixel') } }),
      await shapeRepr('spheres', sphereGrid(3, 3), { props: { clip: clipProps('pixel') } }),
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
