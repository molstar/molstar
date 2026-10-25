/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { background_frag } from './shader/background.frag.js';
import { background_vert } from './shader/background.vert.js';
import { blendBackDpoit_frag } from './shader/blend-back-dpoit.frag.js';
import { blur_frag as bloomBlur_frag } from './shader/bloom/blur.frag.js';
import { composite_frag as bloomComposite_frag } from './shader/bloom/composite.frag.js';
import { luminosity_frag as bloomLuminosity_frag } from './shader/bloom/luminosity.frag.js';
import { cas_frag } from './shader/cas.frag.js';
import { compose_frag as multiSampleCompose_frag } from './shader/compose.frag.js';
import { accumulate_frag as colorSmoothingAccumulate_frag } from './shader/compute/color-smoothing/accumulate.frag.js';
import { accumulate_vert as colorSmoothingAccumulate_vert } from './shader/compute/color-smoothing/accumulate.vert.js';
import { normalize_frag as colorSmoothingNormalize_frag } from './shader/compute/color-smoothing/normalize.frag.js';
import { copy_frag } from './shader/copy.frag.js';
import { cylinders_frag } from './shader/cylinders.frag.js';
import { cylinders_vert } from './shader/cylinders.vert.js';
import { directVolume_frag } from './shader/direct-volume.frag.js';
import { directVolume_vert } from './shader/direct-volume.vert.js';
import { dof_frag } from './shader/dof.frag.js';
import { evaluateDpoit_frag } from './shader/evaluate-dpoit.frag.js';
import { evaluateWboit_frag } from './shader/evaluate-wboit.frag.js';
import { fxaa_frag } from './shader/fxaa.frag.js';
import { gaussianDensity_frag } from './shader/gaussian-density.frag.js';
import { gaussianDensity_vert } from './shader/gaussian-density.vert.js';
import { hiZ_frag } from './shader/hi-z.frag.js';
import { volumeComposite_frag } from './shader/volume-composite.frag.js';
import { reduction_frag } from './shader/histogram-pyramid/reduction.frag.js';
import { sum_frag } from './shader/histogram-pyramid/sum.frag.js';
import { compose_frag as illuminationCompose_frag } from './shader/illumination/compose.frag.js';
import { trace_frag } from './shader/illumination/trace.frag.js';
import { image_frag } from './shader/image.frag.js';
import { image_vert } from './shader/image.vert.js';
import { lines_frag } from './shader/lines.frag.js';
import { lines_vert } from './shader/lines.vert.js';
import { activeVoxels_frag } from './shader/marching-cubes/active-voxels.frag.js';
import { isosurface_frag } from './shader/marching-cubes/isosurface.frag.js';
import { dim_frag } from './shader/marking/dim.frag.js';
import { overlay_frag } from './shader/marking/overlay.frag.js';
import { mesh_frag } from './shader/mesh.frag.js';
import { mesh_vert } from './shader/mesh.vert.js';
import { outlines_frag } from './shader/outlines.frag.js';
import { points_frag } from './shader/points.frag.js';
import { points_vert } from './shader/points.vert.js';
import { postprocessing_frag } from './shader/postprocessing.frag.js';
import { quad_vert } from './shader/quad.vert.js';
import { shadows_frag } from './shader/shadows.frag.js';
import { blend_frag as smaaBlend_frag } from './shader/smaa/blend.frag.js';
import { blend_vert as smaaBlend_vert } from './shader/smaa/blend.vert.js';
import { edges_frag as smaaEdges_frag } from './shader/smaa/edges.frag.js';
import { edges_vert as smaaEdges_vert } from './shader/smaa/edges.vert.js';
import { weights_frag as smaaWeights_frag } from './shader/smaa/weights.frag.js';
import { weights_vert as smaaWeights_vert } from './shader/smaa/weights.vert.js';
import { spheres_frag } from './shader/spheres.frag.js';
import { spheres_vert } from './shader/spheres.vert.js';
import { ssaoBlur_frag } from './shader/ssao-blur.frag.js';
import { ssao_frag } from './shader/ssao.frag.js';
import { text_frag } from './shader/text.frag.js';
import { text_vert } from './shader/text.vert.js';

/** Fragment shader output types, only applicable for WebGL2 */
export type FragOutTypes = { [k in number]: 'vec4' | 'ivec4' };

/** GLSL sources of a program */
export interface GlslProgram {
  readonly vert: string;
  readonly frag: string;
  readonly outTypes: FragOutTypes;
}

const reSingleLineComment = /[ \t]*\/\/.*\n/g;
const reMultiLineComment = /[ \t]*\/\*[\s\S]*?\*\//g;
const reMultipleLinebreaks = /\n{2,}/g;

/** Chunks are composed in TS by the shader modules (see `gl/webgl/shader`), this only strips comments and blank lines */
function normalizeSource(text: string) {
  if (/^\s*#include\b/m.test(text)) {
    throw new Error(`'#include' is not supported, shader sources compose their chunks in TS`);
  }
  return text
    .trim()
    .replace(reSingleLineComment, '\n')
    .replace(reMultiLineComment, '\n')
    .replace(reMultipleLinebreaks, '\n');
}

/**
 * GLSL sources by program name (`ShaderCode.name`), as registered; normalized on first use (`getGlslProgram`): the
 * regexes over every program's sources at module load were most of its import, for programs a session may never draw
 */
const Programs = new Map<string, { source: GlslProgram; normalized?: GlslProgram }>();

/**
 * Registers the GLSL sources of the program `name` (`ShaderCode.name`), replacing earlier ones; for programs
 * generated at runtime (e.g. grid3d) and custom programs
 */
export function registerGlslProgram(name: string, vert: string, frag: string, outTypes: FragOutTypes = {}) {
  Programs.set(name, { source: { vert, frag, outTypes } });
}

/** The GLSL sources of the program `name`, throws for unknown programs */
export function getGlslProgram(name: string): GlslProgram {
  const program = Programs.get(name);
  if (!program) throw new Error(`no GLSL sources for program '${name}'`);
  if (!program.normalized) {
    const { vert, frag, outTypes } = program.source;
    program.normalized = { vert: normalizeSource(vert), frag: normalizeSource(frag), outTypes };
  }
  return program.normalized;
}

export function hasGlslProgram(name: string) {
  return Programs.has(name);
}

/** The names of the registered programs (e.g. for a check that every one compiles) */
export function getGlslProgramNames(): string[] {
  return [...Programs.keys()];
}

// geometry
registerGlslProgram('points', points_vert, points_frag);
registerGlslProgram('spheres', spheres_vert, spheres_frag);
registerGlslProgram('cylinders', cylinders_vert, cylinders_frag);
registerGlslProgram('text', text_vert, text_frag);
registerGlslProgram('lines', lines_vert, lines_frag);
registerGlslProgram('mesh', mesh_vert, mesh_frag);
registerGlslProgram('direct-volume', directVolume_vert, directVolume_frag);
registerGlslProgram('image', image_vert, image_frag);

// passes
registerGlslProgram('background', background_vert, background_frag);
registerGlslProgram('Bloom Luminosity', quad_vert, bloomLuminosity_frag);
registerGlslProgram('Bloom Blur', quad_vert, bloomBlur_frag);
registerGlslProgram('Bloom Composite', quad_vert, bloomComposite_frag);
registerGlslProgram('cas', quad_vert, cas_frag);
registerGlslProgram('dof', quad_vert, dof_frag);
registerGlslProgram('blend-back-dpoit', quad_vert, blendBackDpoit_frag);
registerGlslProgram('evaluate-dpoit', quad_vert, evaluateDpoit_frag);
registerGlslProgram('fxaa', quad_vert, fxaa_frag);
registerGlslProgram('illumination-compose', quad_vert, illuminationCompose_frag);
registerGlslProgram('illumination-multi-sample-compose', quad_vert, multiSampleCompose_frag);
registerGlslProgram('dim', quad_vert, dim_frag);
registerGlslProgram('overlay', quad_vert, overlay_frag);
registerGlslProgram('multi-sample-compose', quad_vert, multiSampleCompose_frag);
registerGlslProgram('outlines', quad_vert, outlines_frag);
registerGlslProgram('postprocessing', quad_vert, postprocessing_frag);
registerGlslProgram('shadows', quad_vert, shadows_frag);
registerGlslProgram('smaa-edges', smaaEdges_vert, smaaEdges_frag);
registerGlslProgram('smaa-weights', smaaWeights_vert, smaaWeights_frag);
registerGlslProgram('smaa-blend', smaaBlend_vert, smaaBlend_frag);
registerGlslProgram('ssao', quad_vert, ssao_frag);
registerGlslProgram('ssao_blur', quad_vert, ssaoBlur_frag);
registerGlslProgram('trace', quad_vert, trace_frag);
registerGlslProgram('evaluate-wboit', quad_vert, evaluateWboit_frag);
registerGlslProgram('copy', quad_vert, copy_frag);
registerGlslProgram('hi-z', quad_vert, hiZ_frag);
registerGlslProgram('volume-composite', quad_vert, volumeComposite_frag);

// compute
registerGlslProgram('gaussian-density', gaussianDensity_vert, gaussianDensity_frag);
registerGlslProgram('color-smoothing-accumulate', colorSmoothingAccumulate_vert, colorSmoothingAccumulate_frag);
registerGlslProgram('normalize', quad_vert, colorSmoothingNormalize_frag);
registerGlslProgram('reduction', quad_vert, reduction_frag, { 0: 'ivec4' });
registerGlslProgram('sum', quad_vert, sum_frag, { 0: 'ivec4' });
registerGlslProgram('active-voxels', quad_vert, activeVoxels_frag);
registerGlslProgram('isosurface', quad_vert, isosurface_frag);
