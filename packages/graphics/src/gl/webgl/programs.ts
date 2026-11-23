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
import { markingDepth_frag } from './shader/marking/depth.frag.js';
import { clipGhost_frag } from './shader/clip-ghost.frag.js';
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
import { gaussianSplats_frag } from './shader/gaussian-splats.frag.js';
import { gaussianSplats_vert } from './shader/gaussian-splats.vert.js';
import { registerGlslProgram } from './program-registry.js';

export {
  type FragOutTypes,
  type GlslProgram,
  getGlslProgram,
  getGlslProgramNames,
  hasGlslProgram,
  registerGlslProgram,
} from './program-registry.js';

// geometry
registerGlslProgram('points', points_vert, points_frag);
registerGlslProgram('spheres', spheres_vert, spheres_frag);
registerGlslProgram('cylinders', cylinders_vert, cylinders_frag);
registerGlslProgram('text', text_vert, text_frag);
registerGlslProgram('lines', lines_vert, lines_frag);
registerGlslProgram('mesh', mesh_vert, mesh_frag);
registerGlslProgram('direct-volume', directVolume_vert, directVolume_frag);
registerGlslProgram('image', image_vert, image_frag);
registerGlslProgram('gaussian-splats', gaussianSplats_vert, gaussianSplats_frag);

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
registerGlslProgram('marking-depth', quad_vert, markingDepth_frag);
registerGlslProgram('clip-ghost', quad_vert, clipGhost_frag);
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
// the first level primed from an occluder grid on WebGPU alone (`HiZPass.setOccluder`): the plain level here
registerGlslProgram('hi-z-occluder', quad_vert, hiZ_frag);
registerGlslProgram('volume-composite', quad_vert, volumeComposite_frag);

// compute
registerGlslProgram('gaussian-density', gaussianDensity_vert, gaussianDensity_frag);
registerGlslProgram('color-smoothing-accumulate', colorSmoothingAccumulate_vert, colorSmoothingAccumulate_frag);
registerGlslProgram('normalize', quad_vert, colorSmoothingNormalize_frag);
registerGlslProgram('reduction', quad_vert, reduction_frag, { 0: 'ivec4' });
registerGlslProgram('sum', quad_vert, sum_frag, { 0: 'ivec4' });
registerGlslProgram('active-voxels', quad_vert, activeVoxels_frag);
registerGlslProgram('isosurface', quad_vert, isosurface_frag);
