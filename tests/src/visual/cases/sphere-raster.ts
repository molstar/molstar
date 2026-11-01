/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { SymmetryOperator } from '@molstar/core/math/geometry/symmetry-operator';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { isDebugMode } from '@molstar/core/util/debug';
import { now } from '@molstar/core/util/now';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { MarkerAction } from '@molstar/core/util/marker-action';
import { ValueCell } from '@molstar/core/util';
import { createClipping } from '@molstar/graphics/geo/geometry/clipping-data';
import { Clipping } from '@molstar/graphics/theme/clipping';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { EveryLoci } from '@molstar/model/model/loci';
import { Structure } from '@molstar/model/model/structure';
import type { PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import type { SceneContext, VisualTestCase } from '../framework/test-case.js';
import { globuleStructure, structureRepr } from '../scenes/structure.js';
import { Impostor, ObliqueView } from './geometry.js';

/**
 * Frames drawn until the sphere raster drew an object (its kernels and the resolve compile in the background, the
 * impostors draw all spheres meanwhile): an error where it never did, on a device with the raster
 */
async function untilRasterDrawn({ graphics, canvas3d }: SceneContext): Promise<string[]> {
  if (!graphics.capabilities.atomicMinMax64 || !graphics.compute) return [];
  const { stats } = graphics;
  // in debug mode the renderer logs why objects are not drawn by the raster: the last of it reported
  const reports: string[] = [];
  const { log } = console;
  console.log = (...args: unknown[]) => {
    const message = args.map(String).join(' ');
    if (message.startsWith('sphere raster')) reports.push(message);
    log.apply(console, args);
  };
  try {
    for (let i = 0; i < 40; ++i) {
      const resolves = stats.sphereRasterResolves;
      await graphics.resources.pending.whenReady();
      canvas3d.requestDraw();
      canvas3d.tick(now(), { isSynchronous: true });
      await graphics.waitForGpuCommandsComplete();
      if (stats.sphereRasterResolves > resolves) return [];
      await new Promise((resolve) => setTimeout(resolve, 25));
    }
  } finally {
    console.log = log;
  }
  const r = canvas3d.props.renderer;
  return [
    `the sphere raster drew no object (${reports.slice(-3).join('; ') || 'nothing logged'}; ${r.sphereRaster} ${r.sphereRasterRadius} ${canvas3d.props.camera.mode} debug ${isDebugMode})`,
  ];
}

/** A grid of `n`³ translations `spacing` apart, centered */
function gridTransforms(n: number, spacing: number) {
  const transforms: Mat4[] = [];
  const o = ((n - 1) * spacing) / 2;
  for (let i = 0; i < n; ++i) {
    for (let j = 0; j < n; ++j) {
      for (let k = 0; k < n; ++k) {
        transforms.push(Mat4.fromTranslation(Mat4(), Vec3.create(i * spacing - o, j * spacing - o, k * spacing - o)));
      }
    }
  }
  return transforms;
}

/** The globule as spacefill, colored per atom */
async function globule(ctx: SceneContext, props: { [k: string]: unknown } = {}) {
  return structureRepr(ctx.graphics, SpacefillRepresentationProvider, await globuleStructure(6000), props as any, {
    name: 'element-index',
  });
}

/**
 * Copies of a small globule on a grid as spacefill, the units merged and culled on the GPU (the raster's input the
 * culler's instances; `merged` false: instanced objects of their own), colored per atom
 */
async function globuleGrid(ctx: SceneContext, props: { [k: string]: unknown } = {}, merged = true) {
  // the copies the structure's own units (`Structure.instances` makes the globule their parent, which the themes read)
  const globule = await globuleStructure(400);
  const builder = Structure.Builder();
  for (const t of gridTransforms(4, 26)) {
    for (const unit of globule.units) builder.addWithOperator(unit, SymmetryOperator.create('copy', t));
  }
  const structure = builder.getStructure();
  return structureRepr(
    ctx.graphics,
    SpacefillRepresentationProvider,
    structure,
    { renderMerged: merged, gpuCulling: merged, ...props } as any,
    { name: 'element-index' },
  );
}

/** A plane through the center, tilted: the instances whose centers are above it clipped whole */
const InstanceClip = {
  variant: 'instance',
  objects: [
    {
      type: 'plane',
      invert: false,
      position: Vec3.create(0, 0, 0),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 30 },
      scale: Vec3.create(1, 1, 1),
      transform: Mat4.identity(),
    },
  ],
};

/** The tilted plane through the center, by pixel: spheres cut, the caps of their solid interiors */
const PixelClip = { ...InstanceClip, variant: 'pixel' };

/** A sphere of 15 Å radius on the surface facing the camera, by pixel: its inside removed (a crater), or its outside (`invert`) */
const SphereClip = (invert: boolean) => ({
  variant: 'pixel',
  objects: [
    {
      type: 'sphere',
      invert,
      position: Vec3.create(9, 11, 22),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 0 },
      scale: Vec3.create(30, 30, 30),
      transform: Mat4.identity(),
    },
  ],
});

/** An ellipsoid on the surface facing the camera, tilted, by pixel: its inside removed */
const EllipsoidClip = {
  variant: 'pixel',
  objects: [
    {
      type: 'sphere',
      invert: false,
      position: Vec3.create(9, 11, 22),
      rotation: { axis: Vec3.create(0, 0, 1), angle: 30 },
      scale: Vec3.create(40, 16, 24),
      transform: Mat4.identity(),
    },
  ],
};

/** An infinite cone from the center, opening toward the camera, by pixel: its inside removed */
const ConeClip = {
  variant: 'pixel',
  objects: [
    {
      type: 'infiniteCone',
      invert: false,
      position: Vec3.create(0, 0, 0),
      rotation: { axis: Vec3.create(1, 0, 0), angle: 200 },
      scale: Vec3.create(1, 1, 1),
      transform: Mat4.identity(),
    },
  ],
};

/**
 * The raster on (`sphereRaster`, its spheres up to `radius` pixels, larger ones impostors), compared with the impostors
 * of the reference (WebGL, which has no raster); devices without the raster draw impostors as well (no `fallback`: it
 * would render the target with the raster off too)
 */
function rasterCase(
  name: string,
  description: string,
  scene: VisualTestCase['scene'],
  options: {
    sphereRaster?: 'on' | 'only';
    radius?: number;
    canvasProps?: PartialCanvas3DProps;
    view?: VisualTestCase['view'];
    interact?: VisualTestCase['interact'];
  } = {},
): VisualTestCase {
  const { sphereRaster = 'on', radius = 6, canvasProps = {} } = options;
  return {
    id: `renderer/sphere-raster/${name}`,
    description,
    requires: Impostor,
    canvasProps: { ...canvasProps, renderer: { ...canvasProps.renderer, sphereRaster, sphereRasterRadius: radius } },
    view: options.view ?? ObliqueView,
    // the canvas's frames, where the raster drew (the image pass's target compiles the resolve again)
    capture: 'canvas',
    scene,
    interact: async (ctx, reprs) => [
      ...(options.interact ? await options.interact(ctx, reprs) : []),
      ...(await untilRasterDrawn(ctx)),
    ],
  };
}

export const SphereRasterCases: VisualTestCase[] = [
  {
    id: 'renderer/sphere-raster/chunks-growing-image',
    description:
      'A sheet of 4096 copies of a globule merged and culled on the GPU, in an image pass: about 192k chunks, more than ' +
      'a merged object starts with, so the first render loses spheres; the read of the counts is pending work, the ' +
      'image rendered again whole.',
    requires: Impostor,
    canvasProps: { renderer: { sphereRaster: 'on', sphereRasterRadius: 64 } },
    // spheres under a pixel: the raster and the impostors pick another nearest sphere at 2.8 % of the pixels, the same
    // with chunks for every instance from the start; lost copies, holes in the sheet, differ at 13 %
    tolerance: { maxDiffRatio: 0.04, maxDiffPixels: 2621 },
    scene: async (ctx) => {
      const globule = await globuleStructure(3000);
      const builder = Structure.Builder();
      for (let i = 0; i < 4096; ++i) {
        const t = Vec3.create(((i % 64) - 31.5) * 26, (Math.floor(i / 64) - 31.5) * 26, 0);
        const operator = SymmetryOperator.create('copy', Mat4.fromTranslation(Mat4(), t));
        for (const unit of globule.units) builder.addWithOperator(unit, operator);
      }
      return [
        await structureRepr(
          ctx.graphics,
          SpacefillRepresentationProvider,
          builder.getStructure(),
          { renderMerged: true, gpuCulling: true } as any,
          { name: 'element-index' },
        ),
      ];
    },
  },
  rasterCase(
    'globule',
    'A globule of 6000 atoms as spacefill, its spheres about 7 px across: those up to 6 px by the raster, the larger ' +
      'near ones as impostors.',
    async (ctx) => [await globule(ctx)],
  ),
  rasterCase(
    'globule-only',
    'The globule with every sphere by the raster (`only`), no impostor drawn.',
    async (ctx) => [await globule(ctx)],
    { sphereRaster: 'only' },
  ),
  rasterCase(
    'globule-near',
    'The globule close up (zoomed): most spheres larger than 6 px as impostors, the far ones by the raster.',
    async (ctx) => [await globule(ctx)],
    { view: { ...ObliqueView, zoom: 2.5 } },
  ),
  rasterCase(
    'merged-grid',
    '64 copies of a globule merged and culled on the GPU: the culler puts the instances the raster draws into regions ' +
      'of their own, the impostors draw the others (the split per instance).',
    async (ctx) => [await globuleGrid(ctx)],
    { radius: 4 },
  ),
  rasterCase(
    'merged-grid-lod',
    'The copies with two levels of detail faded into each other (every 4th sphere, larger, from mid distance): the ' +
      "levels' ranges, scales and fade in the raster as the impostors have them.",
    async (ctx) => [
      await globuleGrid(ctx, {
        lodLevels: [
          { minDistance: 0, maxDistance: 220, overlap: 40, stride: 1, scaleBias: 1 },
          { minDistance: 220, maxDistance: 10000000, overlap: 40, stride: 4, scaleBias: 3 },
        ],
      }),
    ],
    { radius: 4 },
  ),
  rasterCase(
    'clip-instance',
    'The copies as instanced objects (not merged), clipped by instance by a tilted plane: the raster leaves out the ' +
      'instances whose centers are clipped, as the impostors do.',
    async (ctx) => [await globuleGrid(ctx, { clip: InstanceClip }, false)],
    { radius: 4 },
  ),
  rasterCase(
    'merged-clip-instance',
    'The copies merged and culled on the GPU, clipped by instance by a tilted plane, all by the raster (`only`: the ' +
      'clipped copies are the near ones, which `on` leaves to the impostors): the culler leaves out the clipped ' +
      'instances from all regions.',
    async (ctx) => [await globuleGrid(ctx, { clip: InstanceClip })],
    { sphereRaster: 'only' },
  ),
  ...(
    [
      ['clip-pixel', 'the caps of the solid interiors', { clip: PixelClip }],
      [
        'clip-pixel-back',
        'the backs (double sided, no solid interior)',
        { clip: PixelClip, solidInterior: false, doubleSided: true },
      ],
      ['clip-pixel-inverted', 'by an inverted sphere, the caps', { clip: SphereClip(true) }],
      ['clip-pixel-hole', 'holed by a sphere, spheres clipped through, the caps beyond', { clip: SphereClip(false) }],
      ['clip-primitive', 'spheres whole by their centers (`clipPrimitive`)', { clip: PixelClip, clipPrimitive: true }],
      [
        'clip-pixel-ellipsoid',
        'holed by an ellipsoid (the spheres classified by a distance bound)',
        { clip: EllipsoidClip },
      ],
      ['clip-pixel-cone', 'by an infinite cone (the spheres classified by a distance bound)', { clip: ConeClip }],
    ] as const
  ).map(([name, what, props]) =>
    rasterCase(
      name,
      `The globule clipped by pixel, ${what}: per pixel the point a sphere shows as the impostors show it.`,
      async (ctx) => [await globule(ctx, props)],
    ),
  ),
  rasterCase(
    'clip-pixel-layer',
    'The globule clipped by pixel by a tilted plane, every third atom of chain A exempt by a clipping layer (its flags ' +
      'per group): those atoms whole, as the impostors draw them.',
    async (ctx) => {
      const repr = await globule(ctx, { clip: PixelClip });
      const { values } = repr.renderObjects[0];
      const count = values.uGroupCount.ref.value as number;
      createClipping(count, 'groupInstance', values as any);
      const array = values.tClipping.ref.value.array as Uint8Array;
      for (let g = 0; g < count; ++g) array[g] = g % 3 === 0 ? Clipping.Groups.Flag.One : Clipping.Groups.Flag.None;
      ValueCell.update(values.tClipping, values.tClipping.ref.value);
      return [repr];
    },
    // all by the raster: with `on` most of the globule's spheres (about 7 px) are the impostors'
    { sphereRaster: 'only' },
  ),
  rasterCase(
    'merged-clip-pixel',
    'The copies merged and culled on the GPU, clipped by pixel by a tilted plane, all by the raster (`only`): the caps of ' +
      'the instances the plane cuts.',
    async (ctx) => [await globuleGrid(ctx, { clip: PixelClip })],
    { sphereRaster: 'only' },
  ),
  rasterCase(
    'marking-outline',
    'The globule selected (the selection fill and edges drawn over resolved spheres), with outlines.',
    async (ctx) => {
      const repr = await globule(ctx);
      repr.mark(EveryLoci, MarkerAction.Select);
      return [repr];
    },
    {
      canvasProps: {
        postprocessing: {
          outline: { name: 'on', params: PD.getDefaultValues(OutlineParams) },
        },
      },
    },
  ),
  rasterCase('multi-sample', 'The globule multi-sampled: a raster per sample.', async (ctx) => [await globule(ctx)], {
    canvasProps: { multiSample: { mode: 'on', sampleLevel: 2 } },
  }),
  rasterCase(
    'stereo',
    'The globule in side-by-side stereo: the raster draws the left eye (the viewport at the origin), the right eye ' +
      'the impostors.',
    async (ctx) => [await globule(ctx)],
    { canvasProps: { camera: { stereo: { name: 'on', params: { eyeSeparation: 0.062, focus: 10 } } } } },
  ),
];
