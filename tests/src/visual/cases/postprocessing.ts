/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { type PartialCanvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { BloomParams } from '@molstar/graphics/canvas3d/passes/bloom';
import { CasParams } from '@molstar/graphics/canvas3d/passes/cas';
import { DofParams } from '@molstar/graphics/canvas3d/passes/dof';
import { FxaaParams } from '@molstar/graphics/canvas3d/passes/fxaa';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { ShadowParams } from '@molstar/graphics/canvas3d/passes/shadow';
import { SmaaParams } from '@molstar/graphics/canvas3d/passes/smaa';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import { Color } from '@molstar/core/util/color';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { type VisualTestCase } from '../framework/test-case.js';
import { meshBackdrop, meshPrimitives, shapeRepr, sphereGrid } from '../scenes/primitives.js';
import { Impostor, ObliqueView } from './geometry.js';

const ssao = PD.getDefaultValues(SsaoParams);
const shadow = PD.getDefaultValues(ShadowParams);
const outline = PD.getDefaultValues(OutlineParams);
const dof = PD.getDefaultValues(DofParams);
const bloom = PD.getDefaultValues(BloomParams);
const multiScale: (typeof ssao)['multiScale'] = {
  name: 'on',
  params: {
    levels: [
      { radius: 2, bias: 1 },
      { radius: 5, bias: 1 },
      { radius: 8, bias: 1 },
    ],
    nearThreshold: 10,
    farThreshold: 1500,
  },
};

async function primitivesScene() {
  return [
    await shapeRepr('mesh', meshPrimitives()),
    await shapeRepr('backdrop', meshBackdrop(12, -2), { colors: [Color(0xcccccc)] }),
  ];
}

async function transparentScene() {
  return [
    await shapeRepr('mesh', meshPrimitives()),
    await shapeRepr('spheres', sphereGrid(2, 3), { props: { alpha: 0.5 } }),
  ];
}

async function emissiveScene() {
  return [
    await shapeRepr('mesh', meshPrimitives()),
    await shapeRepr('spheres', sphereGrid(2, 3), { size: 0.6, colors: [Color(0xffaa00)], props: { emissive: 0.8 } }),
  ];
}

/** 8x8 PNG, split diagonally into two colors if `color2` is given */
function colorImageUrl(color: string, color2?: string) {
  const canvas = document.createElement('canvas');
  canvas.width = canvas.height = 8;
  const ctx = canvas.getContext('2d')!;
  ctx.fillStyle = color;
  ctx.fillRect(0, 0, 8, 8);
  if (color2) {
    ctx.fillStyle = color2;
    ctx.beginPath();
    ctx.moveTo(8, 0);
    ctx.lineTo(8, 8);
    ctx.lineTo(0, 8);
    ctx.fill();
  }
  return canvas.toDataURL('image/png');
}

function pp(
  id: string,
  postprocessing: PartialCanvas3DProps['postprocessing'],
  scene = primitivesScene,
  extra?: Partial<VisualTestCase>,
): VisualTestCase {
  return {
    id: `postprocessing/${id}`,
    canvasProps: { postprocessing },
    view: ObliqueView,
    scene,
    ...extra,
  };
}

export const PostprocessingCases: VisualTestCase[] = [
  pp('none/primitives', {}),
  pp('occlusion/default', { occlusion: { name: 'on', params: ssao } }, primitivesScene, { tags: ['smoke'] }),
  pp('occlusion/multi-scale', { occlusion: { name: 'on', params: { ...ssao, multiScale } } }),
  // half resolution by default, upsampled with the depth
  pp('occlusion/multi-scale-full-resolution', {
    occlusion: { name: 'on', params: { ...ssao, multiScale, resolutionScale: 1 } },
  }),
  // a depth bias of 0 blurs nothing: the blur passes are skipped
  pp('occlusion/no-blur', { occlusion: { name: 'on', params: { ...ssao, blurDepthBias: 0 } } }),
  pp('occlusion/transparent', { occlusion: { name: 'on', params: ssao } }, transparentScene, { requires: Impostor }),
  pp('shadow/default', { shadow: { name: 'on', params: { ...shadow, steps: 4 } } }),
  pp('outline/default', { outline: { name: 'on', params: outline } }, primitivesScene, { tags: ['smoke'] }),
  pp(
    'outline/transparent',
    { outline: { name: 'on', params: { ...outline, includeTransparent: true } } },
    transparentScene,
    { requires: Impostor },
  ),
  pp('dof/plane', { dof: { name: 'on', params: { ...dof, mode: 'plane' } } }),
  pp('dof/sphere', { dof: { name: 'on', params: { ...dof, mode: 'sphere' } } }),
  // every sample on the pixel: those in focus have no weight
  pp('dof/no-spread', { dof: { name: 'on', params: { ...dof, blurSpread: 0 } } }),
  pp('antialiasing/fxaa', { antialiasing: { name: 'fxaa', params: PD.getDefaultValues(FxaaParams) } }),
  pp('antialiasing/smaa', { antialiasing: { name: 'smaa', params: PD.getDefaultValues(SmaaParams) } }),
  pp('sharpening/cas', { sharpening: { name: 'on', params: PD.getDefaultValues(CasParams) } }),
  pp('bloom/luminosity', { bloom: { name: 'on', params: { ...bloom, mode: 'luminosity', threshold: 0.5 } } }),
  pp('bloom/emissive', { bloom: { name: 'on', params: { ...bloom, mode: 'emissive' } } }, emissiveScene, {
    tags: ['smoke'],
    requires: Impostor,
  }),
  pp('background/horizontal-gradient', {
    background: {
      variant: {
        name: 'horizontalGradient',
        params: { topColor: Color(0x88aadd), bottomColor: Color(0xffffff), ratio: 0.5, coverage: 'viewport' },
      },
    },
  }),
  pp('background/radial-gradient', {
    background: {
      variant: {
        name: 'radialGradient',
        params: { centerColor: Color(0xffffff), edgeColor: Color(0x999999), ratio: 0.5, coverage: 'viewport' },
      },
    },
  }),
  pp('background/skybox', {
    background: {
      variant: {
        name: 'skybox',
        params: {
          faces: {
            name: 'urls',
            params: {
              nx: colorImageUrl('#cc4444'),
              ny: colorImageUrl('#44cc44'),
              nz: colorImageUrl('#4444cc'),
              px: colorImageUrl('#cccc44'),
              py: colorImageUrl('#cc44cc'),
              pz: colorImageUrl('#44cccc'),
            },
          },
          blur: 0,
          rotation: { x: 0, y: 0, z: 0 },
          opacity: 1,
          saturation: 0,
          lightness: 0,
        },
      },
    },
  }),
  pp('background/image', {
    background: {
      variant: {
        name: 'image',
        params: {
          source: { name: 'url', params: colorImageUrl('#88aadd', '#ffffff') },
          blur: 0,
          opacity: 1,
          saturation: 0,
          lightness: 0,
          coverage: 'viewport',
        },
      },
    },
  }),
  pp('background/image-blur', {
    background: {
      variant: {
        name: 'image',
        params: {
          source: { name: 'url', params: colorImageUrl('#88aadd', '#ffffff') },
          blur: 0.3,
          opacity: 1,
          saturation: 0,
          lightness: 0,
          coverage: 'viewport',
        },
      },
    },
  }),
  pp('combined/occlusion-outline-smaa', {
    occlusion: { name: 'on', params: ssao },
    outline: { name: 'on', params: outline },
    antialiasing: { name: 'smaa', params: PD.getDefaultValues(SmaaParams) },
  }),
];
