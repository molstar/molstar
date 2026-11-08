/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { DofParams } from '@molstar/graphics/canvas3d/passes/dof';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import {
  CameraHelperParams,
  type CameraHelperProps,
  isCameraAxesLoci,
} from '@molstar/graphics/canvas3d/helper/camera-helper';
import { Mat4, Vec2, Vec3 } from '@molstar/core/math/linear-algebra';
import type { Clip } from '@molstar/core/util/clip';
import { clipObjectFrame, clipWallsMesh } from '@molstar/graphics/geo/clip-volumes';
import { ModifiersKeys } from '@molstar/core/util/input/input-observer';
import { ScaleLineParams } from '@molstar/graphics/canvas3d/helper/scale-line';
import { LinesBuilder } from '@molstar/graphics/geo/geometry/lines/lines-builder';
import { PointsBuilder } from '@molstar/graphics/geo/geometry/points/points-builder';
import { Representation } from '@molstar/graphics/repr/representation';
import { transformHandles } from '@molstar/graphics/util/gizmo/transform';
import { cartesian, mergeRequires, variants } from '../framework/matrix.js';
import type { SceneContext, VisualTestCase } from '../framework/test-case.js';
import { meshBox, meshTorus, shapeRepr, sphereGrid } from '../scenes/primitives.js';
import { Impostor, ObliqueView } from './geometry.js';

/** How the scene's depth applies to the overlay (`occludedAlpha`) */
const occlusions = variants<{ occludedAlpha: number }>({
  'on-top': { occludedAlpha: 1 },
  'x-ray': { occludedAlpha: 0.25 },
  'depth-tested': { occludedAlpha: 0 },
});

/** A gray box in the scene, a torus (mesh) and spheres (impostors) through it in the overlay */
function overlayScene(props: { occludedAlpha: number }, spheres: boolean) {
  return async (): Promise<Representation.Any[]> => {
    const box = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] });
    const torus = await shapeRepr('torus', meshTorus(3), { colors: [Color(0x377eb8)], props });
    torus.setState({ layer: 'overlay' });
    if (!spheres) return [box, torus];
    const balls = await shapeRepr('spheres', sphereGrid(2, 4), { colors: [Color(0xe41a1c)], props });
    balls.setState({ layer: 'overlay' });
    return [box, torus, balls];
  };
}

const clipObject = (
  type: Clip.Props['objects'][number]['type'],
  position: Vec3,
  scale: Vec3,
  rotation = { axis: Vec3.create(1, 0, 0), angle: 0 },
): Clip.Props['objects'][number] => ({ type, invert: false, position, rotation, scale, transform: Mat4.identity() });

/** The clip objects whose walls lie on the caps they cut into a solid box */
const clipWallObjects = variants<Clip.Props['objects'][number]>({
  plane: clipObject('plane', Vec3(), Vec3.create(1, 1, 1), { axis: Vec3.create(0, 0, 1), angle: 30 }),
  sphere: clipObject('sphere', Vec3.create(2, 2, 2), Vec3.create(4, 4, 4)),
  cube: clipObject('cube', Vec3.create(2, 2, 2), Vec3.create(4, 4, 4), { axis: Vec3.create(0, 1, 0), angle: 25 }),
  cylinder: clipObject('cylinder', Vec3.create(2, 2, 2), Vec3.create(4, 4, 4)),
  'infinite-cone': clipObject('infiniteCone', Vec3.create(0, 0, 6), Vec3.create(1, 1, 1)),
});

/** A solid box cut by `object`, and the walls on its cap: in the scene, or in the overlay as a gizmo's body */
function clipWallsScene(object: Clip.Props['objects'][number], overlay: boolean) {
  return async (): Promise<Representation.Any[]> => {
    const box = await shapeRepr('box', meshBox(6), {
      colors: [Color(0x999999)],
      props: { solidInterior: true, clip: { variant: 'pixel', objects: [object] } },
    });
    const walls = await shapeRepr('walls', clipWallsMesh(object, 6), {
      colors: [Color(0x377eb8)],
      props: { alpha: 0.4, ignoreLight: true, ...(overlay ? { occludedAlpha: 0.25 } : { depthOffset: 0.02 }) },
    });
    walls.setState({
      transform: clipObjectFrame(object, { center: Vec3(), radius: 6 } as Parameters<typeof clipObjectFrame>[1]),
      ...(overlay ? { layer: 'overlay' as const } : {}),
    });
    return [box, walls];
  };
}

/** A gray box in the scene, a screen-sized torus (mesh) through it in the overlay, x-ray */
async function screenSizeScene(): Promise<Representation.Any[]> {
  const box = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] });
  const torus = await shapeRepr('torus', meshTorus(1), {
    colors: [Color(0x377eb8)],
    props: { occludedAlpha: 0.25, screenSize: 20 },
  });
  torus.setState({ layer: 'overlay' });
  return [box, torus];
}

/** A gray box in the scene, a circle in the overlay facing the camera (a `ring` point): round from any side */
async function facingScene(): Promise<Representation.Any[]> {
  const box = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] });
  const center = PointsBuilder.create(1, 1);
  center.add(0, 0, 0, 0);
  const circle = await shapeRepr('circle', center.getPoints(), {
    colors: [Color(0x377eb8)],
    size: 120,
    props: {
      occludedAlpha: 0.25,
      sizeFactor: 1,
      pointSizeAttenuation: false,
      pointStyle: 'ring' as const,
      pointRingWidth: 0.1,
    },
  });
  circle.setState({ layer: 'overlay' });
  return [box, circle];
}

/**
 * The transform handles' parts at rest over the gray box, as the controller draws them at `HandlesScreenSize` CSS pixels
 * a unit: the heads, squares and front rings meshes, the center and roll circles `ring` points, no hit areas
 */
const HandlesScreenSize = 60;
async function handlePartsScene(): Promise<Representation.Any[]> {
  const box = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] });
  const { parts, groups } = transformHandles(1);
  const colors = groups.map((g) => g.color);
  const out: Representation.Any[] = [box];
  for (const part of parts) {
    if (part.role !== 'handle' || part.when) continue;
    const g = part.geometry;
    const p = part.point;
    const repr = p
      ? await shapeRepr('part', g, {
          colors,
          size: 2 * p.radius * HandlesScreenSize,
          props: {
            emissive: 0.35,
            sizeFactor: 1,
            pointSizeAttenuation: false,
            pointStyle: p.style ?? 'circle',
            pointRingWidth: p.ringWidth ?? 1,
            pointDepthOffset: (p.depth ?? 0) * HandlesScreenSize,
          },
        })
      : await shapeRepr('part', g, {
          colors,
          props: { ignoreLight: true, emissive: 0.35, screenSize: HandlesScreenSize },
        });
    repr.setState({ layer: 'overlay' });
    out.push(repr);
  }
  return out;
}

/** The transform handles inside the gray box, on top, as gizmos draw them: unlit, glowing */
function handlesScene() {
  return async (): Promise<Representation.Any[]> => {
    const box = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] });
    const { handles, groups } = transformHandles(4);
    const gizmo = await shapeRepr('handles', handles, {
      colors: groups.map((g) => g.color),
      props: { ignoreLight: true, emissive: 0.35 },
    });
    gizmo.setState({ layer: 'overlay' });
    return [box, gizmo];
  };
}

/** The pointer over the corner axes (bottom-left): a part of them picked near their origin, hovered as Canvas3D does */
async function hoverAxes({ graphics, canvas3d }: SceneContext): Promise<string[]> {
  const height = canvas3d.camera.viewport.height / graphics.pixelRatio;
  // their origin inset by their bounds' radius (41 CSS pixels at the default scale)
  for (let dy = 0; dy <= 24; dy += 4) {
    for (let dx = 0; dx <= 24; dx += 4) {
      const data = await canvas3d.asyncIdentify(Vec2.create(41 + dx, height - 41 - dy))?.wait();
      const current = data && canvas3d.getLoci(data.id);
      if (!current || !isCameraAxesLoci(current.loci)) continue;
      canvas3d.interaction.hover.next({ current, buttons: 0, button: 0, modifiers: ModifiersKeys.None });
      return [];
    }
  }
  return ['no part of the corner axes picked near their origin'];
}

/**
 * A clip gizmo's shape over the gray box, as the clipping gizmos draw it: a 6 Å box's edges as lines and its corners
 * as round points, a set size on screen
 */
async function clipShapeScene(): Promise<Representation.Any[]> {
  const box = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] });
  const edges = LinesBuilder.create(16, 16);
  const corners = PointsBuilder.create(8, 8);
  const h = 3;
  for (const a of [0, 1, 2]) {
    const [b, c] = [0, 1, 2].filter((x) => x !== a);
    for (const sb of [-1, 1]) {
      for (const sc of [-1, 1]) {
        const from = [0, 0, 0];
        const to = [0, 0, 0];
        from[a] = -h;
        to[a] = h;
        from[b] = to[b] = sb * h;
        from[c] = to[c] = sc * h;
        edges.add(from[0], from[1], from[2], to[0], to[1], to[2], 0);
      }
    }
  }
  for (const x of [-1, 1]) for (const y of [-1, 1]) for (const z of [-1, 1]) corners.add(x * h, y * h, z * h, 0);
  const lines = await shapeRepr('edges', edges.getLines(), {
    colors: [Color(0x9db8e0)],
    size: 2.5,
    props: { sizeFactor: 1, lineSizeAttenuation: false },
  });
  const points = await shapeRepr('corners', corners.getPoints(), {
    colors: [Color(0x9db8e0)],
    size: 12,
    props: { sizeFactor: 1, pointSizeAttenuation: false, pointStyle: 'circle' },
  });
  for (const r of [lines, points]) r.setState({ layer: 'overlay' });
  return [box, lines, points];
}

/** The corner axes on at `location`; labels and an offset in the top-right */
function axesProps(location: (typeof Representation.ScreenEdgeAnchors)[number]): CameraHelperProps {
  const params = PD.getDefaultValues((CameraHelperParams.axes.map('on') as PD.Group<PD.Params>).params);
  const labels = location === 'top-right';
  return {
    axes: {
      name: 'on',
      params: {
        ...(params as CameraHelperProps['axes']['params']),
        location,
        showLabels: labels,
        locationOffsetX: labels ? 12 : 0,
        locationOffsetY: labels ? 12 : 0,
      },
    },
  } as CameraHelperProps;
}

/** DOF that blurs all but a thin slice at the scene's center */
const StrongDof = { name: 'on', params: { ...PD.getDefaultValues(DofParams), blurSize: 16, PPM: 1 } } as const;

const passes = variants<Pick<VisualTestCase, 'contextProps' | 'canvasProps' | 'capture'>>(
  {
    blended: {},
    wboit: { contextProps: { transparency: 'wboit' } },
    dpoit: { contextProps: { transparency: 'dpoit' } },
    postprocessing: {
      canvasProps: {
        postprocessing: {
          occlusion: { name: 'on', params: PD.getDefaultValues(SsaoParams) },
          outline: { name: 'on', params: PD.getDefaultValues(OutlineParams) },
        },
      },
    },
    'multi-sample': { canvasProps: { multiSample: { mode: 'on', sampleLevel: 2 } } },
    illumination: { canvasProps: { illumination: { enabled: true, maxIterations: 3 } } },
    canvas: { capture: 'canvas' },
    // the scene blurred, the overlay after DOF: sharp
    dof: { canvasProps: { postprocessing: { dof: StrongDof } } },
    'dof-illumination': {
      canvasProps: { postprocessing: { dof: StrongDof }, illumination: { enabled: true, maxIterations: 3 } },
    },
    'dof-canvas': { capture: 'canvas', canvasProps: { postprocessing: { dof: StrongDof } } },
  },
  {
    wboit: { requires: ['wboit'] },
    dpoit: { requires: ['dpoit'] },
    illumination: { requires: ['illumination'] },
    'dof-illumination': { requires: ['illumination'] },
  },
);

export const OverlayCases: VisualTestCase[] = [
  ...occlusions.map(
    (occlusion): VisualTestCase => ({
      id: `overlay/occlusion/${occlusion.name}`,
      description: 'The overlay through a box of the scene: on top, x-ray, depth tested, x-ray without color.',
      tags: occlusion.name === 'x-ray' ? ['smoke'] : [],
      requires: Impostor,
      view: ObliqueView,
      scene: overlayScene(occlusion.value, true),
    }),
  ),
  ...cartesian(
    passes.filter((p) => p.name !== 'blended'),
    occlusions.filter((o) => o.name === 'x-ray'),
    (pass, occlusion): VisualTestCase => ({
      id: `overlay/passes/${pass.name}`,
      description: 'The x-ray overlay drawn by each path: over the finished image, before antialiasing.',
      requires: mergeRequires(pass.requires, Impostor),
      ...pass.value,
      view: ObliqueView,
      scene: overlayScene(occlusion.value, true),
    }),
  ),
  ...([1, 0.5] as const).map(
    (zoom): VisualTestCase => ({
      id: `overlay/screen-size/zoom-${zoom}`,
      description: 'A screen-sized torus (20 pixels per unit) over the box: the same size at both zooms, the box not.',
      view: { ...ObliqueView, zoom },
      scene: screenSizeScene,
    }),
  ),
  {
    id: 'overlay/screen-size/stereo',
    description: 'The screen-sized torus in each eye of side-by-side stereo, by the eye projection.',
    capture: 'canvas',
    size: { width: 384, height: 192 },
    canvasProps: { camera: { stereo: { name: 'on', params: { eyeSeparation: 0.062, focus: 10 } } } },
    view: ObliqueView,
    scene: screenSizeScene,
  },
  {
    id: 'overlay/handles',
    description: 'The transform handles unlit, glowing, on top of the box.',
    tags: ['smoke'],
    view: ObliqueView,
    scene: handlesScene(),
  },
  {
    id: 'overlay/handles/parts',
    description:
      'The transform handles at rest, a set size on screen: heads, squares and front rings meshes, the center and roll circles ring points.',
    view: ObliqueView,
    scene: handlePartsScene,
  },
  ...([1, 0.5] as const).map(
    (zoom): VisualTestCase => ({
      id: `overlay/clip-shape/zoom-${zoom}`,
      description: 'A clip gizmo shape: edges as lines, corners as round points (a set size).',
      view: { ...ObliqueView, zoom },
      scene: clipShapeScene,
    }),
  ),
  {
    id: 'overlay/facing/side',
    description: 'A circle facing the camera (a ring point) seen from the side: round, not a line.',
    view: { direction: [-1, 0, 0] },
    scene: facingScene,
  },
  {
    id: 'overlay/facing/stereo',
    description: 'The facing circle in each eye of side-by-side stereo: round in each.',
    capture: 'canvas',
    size: { width: 384, height: 192 },
    canvasProps: { camera: { stereo: { name: 'on', params: { eyeSeparation: 0.062, focus: 10 } } } },
    view: { direction: [-1, 0, 0] },
    scene: facingScene,
  },
  ...(['bottom-left', 'bottom-right', 'top-left', 'top-right'] as const).map(
    (location): VisualTestCase => ({
      id: `screen/axes/${location}`,
      description: 'The corner axes in the screen layer, over the box; labels and an offset in the top-right.',
      canvasProps: { camera: { helper: axesProps(location) } },
      view: ObliqueView,
      scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
    }),
  ),
  {
    id: 'screen/axes/canvas',
    description: 'The corner axes captured from the canvas (Canvas3D frame) instead of an image pass.',
    capture: 'canvas',
    canvasProps: { camera: { helper: axesProps('bottom-left') } },
    view: ObliqueView,
    scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
  },
  {
    id: 'screen/axes/hovered',
    description: 'The corner axes with the pointer over them: at full alpha (else at their dimmed alpha).',
    capture: 'canvas',
    canvasProps: { camera: { helper: axesProps('bottom-left') } },
    view: ObliqueView,
    scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
    interact: hoverAxes,
  },
  {
    id: 'screen/scale-line/over-the-scene',
    description: 'The scale line at the bottom middle, its grey label over the grey box: readable by its outline.',
    canvasProps: {
      camera: {
        mode: 'orthographic',
        helper: {
          axes: { name: 'off', params: {} },
          scale: { name: 'on', params: { ...PD.getDefaultValues(ScaleLineParams), location: 'bottom-middle' } },
        },
      },
    },
    view: ObliqueView,
    scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
  },
  {
    id: 'screen/middle',
    description: 'The axes at the top middle, the scale line at the bottom middle (centered on it).',
    canvasProps: {
      camera: {
        mode: 'orthographic',
        helper: {
          ...axesProps('top-middle'),
          scale: { name: 'on', params: { ...PD.getDefaultValues(ScaleLineParams), location: 'bottom-middle' } },
        },
      },
    },
    // zoomed out: the grey label not over the grey box (no contrast there)
    view: { ...ObliqueView, zoom: 0.6 },
    scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
  },
  ...(['image', 'canvas'] as const).map(
    (capture): VisualTestCase => ({
      id: `screen/scale-line/${capture}`,
      description: 'The scale line (bottom-middle) over the box, orthographic: a 1, 2, 5 × 10ⁿ Å length and its label.',
      capture: capture === 'image' ? 'image-pass' : 'canvas',
      canvasProps: {
        camera: {
          mode: 'orthographic',
          helper: {
            axes: { name: 'off', params: {} },
            scale: { name: 'on', params: PD.getDefaultValues(ScaleLineParams) },
          },
        },
      },
      view: ObliqueView,
      scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
    }),
  ),
  {
    id: 'screen/axes/stereo',
    description: 'The corner axes in each eye of side-by-side stereo.',
    capture: 'canvas',
    size: { width: 384, height: 192 },
    canvasProps: {
      camera: {
        helper: axesProps('bottom-left'),
        stereo: { name: 'on', params: { eyeSeparation: 0.062, focus: 10 } },
      },
    },
    view: ObliqueView,
    scene: async () => [await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] })],
  },
  ...cartesian(
    variants<{ alpha: number }>({ picked: { alpha: 0.8 }, faint: { alpha: 0.2 } }),
    // where wboit or dpoit is not supported, the target draws blended: so does the reference
    variants<Pick<VisualTestCase, 'contextProps' | 'fallback'>>({
      blended: {},
      wboit: {
        contextProps: { transparency: 'wboit' },
        fallback: { requires: ['wboit'], overrides: { contextProps: { transparency: 'blended' } } },
      },
      dpoit: {
        contextProps: { transparency: 'dpoit' },
        fallback: { requires: ['dpoit'], overrides: { contextProps: { transparency: 'blended' } } },
      },
    }),
    (box, pass) => ({
      id: `overlay/occlusion/transparent/${box.name}/${pass.name}`,
      description:
        'A depth-tested overlay torus through a transparent box: hidden by a box opaque enough to be picked, shown through a fainter one.',
      ...pass.value,
      view: ObliqueView,
      scene: async () => {
        const scene = await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)], props: box.value });
        const torus = await shapeRepr('torus', meshTorus(3), {
          colors: [Color(0x377eb8)],
          props: { occludedAlpha: 0 },
        });
        torus.setState({ layer: 'overlay' });
        return [scene, torus];
      },
    }),
  ),
  ...cartesian(clipWallObjects, variants<boolean>({ scene: false, overlay: true }), (object, layer) => ({
    id: `overlay/clip-walls/${object.name}/${layer.name}`,
    description: 'The walls of a clip object on the solid-interior cap it cuts: even, no z-fighting.',
    requires: Impostor,
    view: ObliqueView,
    scene: clipWallsScene(object.value, layer.value),
  })),
  {
    id: 'overlay/occlusion/mesh-only',
    description: 'The x-ray overlay without impostors, on every profile (no frag depth needed).',
    view: ObliqueView,
    scene: overlayScene(occlusions.find((o) => o.name === 'x-ray')!.value, false),
  },
];
