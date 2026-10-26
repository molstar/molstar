/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { Canvas3DProps } from '@molstar/graphics/canvas3d/canvas3d';
import { BloomParams } from '@molstar/graphics/canvas3d/passes/bloom';
import { CasParams } from '@molstar/graphics/canvas3d/passes/cas';
import { DofParams } from '@molstar/graphics/canvas3d/passes/dof';
import { FxaaParams } from '@molstar/graphics/canvas3d/passes/fxaa';
import { IlluminationParams } from '@molstar/graphics/canvas3d/passes/illumination';
import { OutlineParams } from '@molstar/graphics/canvas3d/passes/outline';
import { PostprocessingParams, type PostprocessingProps } from '@molstar/graphics/canvas3d/passes/postprocessing';
import { ShadowParams } from '@molstar/graphics/canvas3d/passes/shadow';
import { SmaaParams } from '@molstar/graphics/canvas3d/passes/smaa';
import { SsaoParams } from '@molstar/graphics/canvas3d/passes/ssao';
import type { RepresentationProvider } from '@molstar/graphics/repr/representation';
import { BallAndStickRepresentationProvider } from '@molstar/graphics/repr/structure/representation/ball-and-stick';
import { CartoonRepresentationProvider } from '@molstar/graphics/repr/structure/representation/cartoon';
import { GaussianSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/gaussian-surface';
import { MolecularSurfaceRepresentationProvider } from '@molstar/graphics/repr/structure/representation/molecular-surface';
import { SpacefillRepresentationProvider } from '@molstar/graphics/repr/structure/representation/spacefill';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { structureRepr } from '../visual/scenes/structure.js';
import type { BenchmarkCase } from './framework.js';
import { getStructure, LargeSize, type RenderInput, renderCase } from './render.js';

// GPU time per pass (`--gpu-timing`), not in the default runs (`--filter 'passes/*'`): each pass on its own over a
// scene filling a 2560x1440 canvas, so that full-screen passes stand out of the noise; the geometry cases draw the
// opaque pass alone (no postprocessing, no antialiasing), the others add one pass to `passes/spheres`

const ssao = PD.getDefaultValues(SsaoParams);

/** postprocessing with every pass off, `pass` on */
function postprocessing(pass: Partial<PostprocessingProps> = {}): Partial<Canvas3DProps> {
  return {
    postprocessing: {
      ...PD.getDefaultValues(PostprocessingParams),
      occlusion: { name: 'off', params: {} },
      antialiasing: { name: 'off', params: {} },
      bloom: { name: 'off', params: {} },
      ...pass,
    },
  };
}

/** 7qpd zoomed in to fill the canvas, as `providers` */
function scene(
  providers: { provider: RepresentationProvider<any, any, any>; props?: object }[],
  input: Partial<RenderInput> = {},
): () => Promise<RenderInput> {
  return async () => {
    const s = await getStructure();
    return {
      size: LargeSize,
      zoom: 0.6,
      scene: (g) => Promise.all(providers.map(({ provider, props }) => structureRepr(g, provider, s, props ?? {}))),
      props: postprocessing(),
      ...input,
    };
  };
}

const spheres = [{ provider: SpacefillRepresentationProvider }];
/** a cartoon inside a transparent molecular surface: the transparent passes over opaque geometry */
const transparent = [
  { provider: CartoonRepresentationProvider },
  { provider: MolecularSurfaceRepresentationProvider, props: { alpha: 0.5 } },
];

const pass = (id: string, description: string, props: Partial<PostprocessingProps>) =>
  renderCase(`passes/${id}`, `spacefill with ${description}`, scene(spheres, { props: postprocessing(props) }));

export const PassCases: BenchmarkCase<any>[] = [
  renderCase('passes/mesh', 'cartoon (meshes)', scene([{ provider: CartoonRepresentationProvider }])),
  renderCase('passes/spheres', 'spacefill (sphere impostors)', scene(spheres)),
  // the impostor's parts apart: its shading, and its ray-sphere hit and depth write
  renderCase(
    'passes/spheres-ignore-light',
    'spacefill (sphere impostors) without lighting',
    scene([{ provider: SpacefillRepresentationProvider, props: { ignoreLight: true } }]),
  ),
  renderCase(
    'passes/spheres-approximate',
    'spacefill as approximate impostors (no ray-sphere hit, no depth write)',
    scene([{ provider: SpacefillRepresentationProvider, props: { approximate: true } }]),
  ),
  renderCase(
    'passes/cylinders',
    'ball-and-stick (sphere and cylinder impostors)',
    scene([{ provider: BallAndStickRepresentationProvider }]),
  ),
  renderCase(
    'passes/gaussian-surface',
    'Gaussian surface (GPU texture meshes)',
    scene([{ provider: GaussianSurfaceRepresentationProvider }]),
  ),
  renderCase('passes/marking', 'spacefill, everything highlighted', scene(spheres, { highlight: true })),
  pass('ssao', 'SSAO', { occlusion: { name: 'on', params: ssao } }),
  pass('ssao-multiscale', 'multi-scale SSAO (four levels)', {
    occlusion: {
      name: 'on',
      params: { ...ssao, multiScale: { name: 'on', params: SsaoParams.multiScale.map('on').defaultValue } },
    },
  }),
  pass('outline', 'outlines', { outline: { name: 'on', params: PD.getDefaultValues(OutlineParams) } }),
  pass('shadow', 'shadows (4 steps)', {
    shadow: { name: 'on', params: { ...PD.getDefaultValues(ShadowParams), steps: 4 } },
  }),
  pass('dof', 'depth of field', { dof: { name: 'on', params: PD.getDefaultValues(DofParams) } }),
  pass('bloom', 'bloom (luminosity)', {
    bloom: { name: 'on', params: { ...PD.getDefaultValues(BloomParams), mode: 'luminosity', threshold: 0.5 } },
  }),
  pass('fxaa', 'FXAA', { antialiasing: { name: 'fxaa', params: PD.getDefaultValues(FxaaParams) } }),
  pass('smaa', 'SMAA', { antialiasing: { name: 'smaa', params: PD.getDefaultValues(SmaaParams) } }),
  pass('cas', 'sharpening (CAS)', { sharpening: { name: 'on', params: PD.getDefaultValues(CasParams) } }),
  ...(['blended', 'wboit', 'dpoit'] as const).map((mode) =>
    renderCase(
      `passes/transparency-${mode}`,
      `cartoon in a transparent molecular surface, ${mode} transparency`,
      scene(transparent, { transparency: mode }),
    ),
  ),
  renderCase(
    'passes/illumination',
    'cartoon with illumination: the first tracing iteration per frame (fixed quality)',
    scene([{ provider: CartoonRepresentationProvider }], {
      props: {
        ...postprocessing(),
        illumination: { ...PD.getDefaultValues(IlluminationParams), enabled: true, rendersPerFrame: [4, 4] },
      },
    }),
  ),
];
