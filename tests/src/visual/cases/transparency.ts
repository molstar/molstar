/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Canvas3DContext } from '@molstar/graphics/canvas3d/canvas3d';
import { Mat4, Vec3 } from '@molstar/core/math/linear-algebra';
import { Color } from '@molstar/core/util/color';
import { cartesian, mergeRequires, variants } from '../framework/matrix.js';
import { type SceneContext, type VisualTestCase } from '../framework/test-case.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { meshBox, meshTorus, shapeRepr, sphereGrid } from '../scenes/primitives.js';
import { blobVolume, volumeRepr } from '../scenes/volume.js';
import { DirectVolumeRepresentationProvider } from '@molstar/graphics/repr/volume/direct-volume';
import { Impostor, ObliqueView } from './geometry.js';
import { Representation } from '@molstar/graphics/repr/representation';

type Transparency = Canvas3DContext.Props['transparency'];

const modes = variants<{ transparency: Transparency; dpoitIterations?: number }>(
  {
    blended: { transparency: 'blended' },
    wboit: { transparency: 'wboit' },
    dpoit: { transparency: 'dpoit' },
    'dpoit-4': { transparency: 'dpoit', dpoitIterations: 4 },
  },
  {
    wboit: { requires: ['wboit'] },
    dpoit: { requires: ['dpoit'] },
    'dpoit-4': { requires: ['dpoit'] },
  },
);

/**
 * A direct volume behind a transparent box (between it and the camera of \`ObliqueView\`): the volumes are drawn before
 * the transparent objects, which cover them
 */
async function volumeBehindScene({ graphics }: SceneContext) {
  return [
    await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(32), {}),
    await shapeRepr('box', meshBox(10), {
      colors: [Color(0x377eb8)],
      transforms: [Mat4.fromTranslation(Mat4(), Vec3.create(26, 20, 28))],
      props: { alpha: 0.5 },
    }),
  ];
}

const scenes = variants<(ctx: SceneContext) => Promise<Representation.Any[]>>(
  {
    overlap: async () => [
      await shapeRepr('box', meshBox(5), { colors: [Color(0x999999)] }),
      await shapeRepr('spheres', sphereGrid(2, 3), { props: { alpha: 0.5 } }),
      await shapeRepr('torus', meshTorus(3), { colors: [Color(0x377eb8)], props: { alpha: 0.6 } }),
    ],
    'low-alpha': async () => [await shapeRepr('spheres', sphereGrid(3, 2.5), { props: { alpha: 0.25 } })],
    backfaces: async () => [
      await shapeRepr('torus', meshTorus(3), { props: { alpha: 0.5, transparentBackfaces: 'on', doubleSided: true } }),
    ],
    /**
     * a direct volume around an opaque box: each mode draws the volume its own way (at half resolution weighted blended
     * targets of its own or the volume pass), cut at the box
     */
    volume: async ({ graphics }) => [
      await volumeRepr(graphics, DirectVolumeRepresentationProvider, blobVolume(32), {}),
      await shapeRepr('box', meshBox(8), {
        colors: [Color(0x999999)],
        transforms: [Mat4.fromTranslation(Mat4(), Vec3.create(16, 16, 16))],
      }),
    ],
    'volume-behind': volumeBehindScene,
    /** caps of the clipped torus, stencil-marked and filled in the transparent passes */
    'solid-interior': async () => [
      await shapeRepr('torus', meshTorus(3), {
        props: {
          alpha: 0.5,
          transparentBackfaces: 'on',
          doubleSided: true,
          solidInterior: true,
          clip: {
            variant: 'pixel',
            objects: [
              {
                type: 'plane',
                invert: false,
                position: Vec3.create(0, 0, 0),
                rotation: { axis: Vec3.create(0, 1, 0), angle: 60 },
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
    overlap: { requires: Impostor },
    'low-alpha': { requires: Impostor },
    // caps need frag depth
    'solid-interior': { requires: Impostor },
  },
);

export const TransparencyCases: VisualTestCase[] = [
  ...cartesian(modes, scenes, (mode, scene) => ({
    id: `transparency/${mode.name}/${scene.name}`,
    tags: scene.name === 'overlap' ? ['smoke'] : [],
    requires: mergeRequires(mode.requires, scene.requires),
    contextProps: { transparency: mode.value.transparency },
    canvasProps: mode.value.dpoitIterations ? { dpoitIterations: mode.value.dpoitIterations } : undefined,
    view: ObliqueView,
    scene: scene.value,
  })),
  ...(['blended', 'wboit', 'dpoit'] as const).flatMap((transparency): VisualTestCase[] => [
    {
      id: `transparency/${transparency}/volume-behind-full-resolution`,
      description: 'The volume behind the transparent box at the full resolution, as the transparency mode draws it.',
      requires: transparency === 'blended' ? undefined : [transparency],
      contextProps: { transparency },
      canvasProps: { volumeResolutionScale: 1 },
      view: ObliqueView,
      scene: volumeBehindScene,
    },
    {
      id: `transparency/${transparency}/volume-behind-postprocessing`,
      description:
        'The volume behind the transparent box, drawn into the transparent layer that postprocessing blends.',
      requires: transparency === 'blended' ? undefined : [transparency],
      contextProps: { transparency },
      canvasProps: { postprocessing: { occlusion: { name: 'on', params: PD.getDefaultValues(SsaoParams) } } },
      view: ObliqueView,
      scene: volumeBehindScene,
    },
  ]),
  {
    id: 'transparency/illumination/volume-behind',
    description: 'The volume behind the transparent box with the traced shading.',
    requires: ['illumination'],
    canvasProps: { illumination: { enabled: true, maxIterations: 3 } },
    view: ObliqueView,
    scene: volumeBehindScene,
  },
  {
    id: 'transparency/background/transparent',
    description: 'Alpha channel is compared in memory; PNG snapshots of non-opaque images are lossy.',
    requires: Impostor,
    canvasProps: { transparentBackground: true },
    view: ObliqueView,
    scene: async () => [await shapeRepr('spheres', sphereGrid(2, 3), { props: { alpha: 0.5 } })],
  },
  ...[false, true].map(
    (fog): VisualTestCase => ({
      id: `transparency/background/transparent-dpoit${fog ? '-fog' : ''}`,
      description: 'Dual depth peeling keeps opaque and transparent colors apart until the blend, also with fog.',
      requires: mergeRequires(['dpoit'], Impostor),
      contextProps: { transparency: 'dpoit' },
      canvasProps: {
        transparentBackground: true,
        ...(fog ? { cameraFog: { name: 'on', params: { intensity: 50 } } } : {}),
      },
      view: ObliqueView,
      scene: scenes.find((s) => s.name === 'overlap')!.value,
    }),
  ),
];
