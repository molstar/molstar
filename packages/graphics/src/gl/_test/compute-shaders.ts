/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import * as background from '../../canvas3d/passes/background.js';
import * as bloom from '../../canvas3d/passes/bloom.js';
import * as cas from '../../canvas3d/passes/cas.js';
import * as dof from '../../canvas3d/passes/dof.js';
import * as dpoit from '../../canvas3d/passes/dpoit.js';
import * as fxaa from '../../canvas3d/passes/fxaa.js';
import * as illumination from '../../canvas3d/passes/illumination.js';
import * as marking from '../../canvas3d/passes/marking.js';
import * as multiSample from '../../canvas3d/passes/multi-sample.js';
import * as outline from '../../canvas3d/passes/outline.js';
import * as postprocessing from '../../canvas3d/passes/postprocessing.js';
import * as shadow from '../../canvas3d/passes/shadow.js';
import * as smaa from '../../canvas3d/passes/smaa.js';
import * as ssao from '../../canvas3d/passes/ssao.js';
import * as tracing from '../../canvas3d/passes/tracing.js';
import * as wboit from '../../canvas3d/passes/wboit.js';
import * as volume from '../../canvas3d/passes/volume.js';
import * as gaussianDensity from '../../geo/gaussian-density/gpu.js';
import * as colorSmoothing from '../../geo/geometry/texture-mesh/color-smoothing.js';
import * as histogramReduction from '../webgl/compute/histogram-pyramid/reduction.js';
import * as histogramSum from '../webgl/compute/histogram-pyramid/sum.js';
import * as hiZ from '../compute/hi-z.js';
import * as activeVoxels from '../webgl/compute/marching-cubes/active-voxels.js';
import * as isosurface from '../webgl/compute/marching-cubes/isosurface.js';
import * as computeUtil from '../compute/util.js';
import type { RenderableSchema } from '../renderable/schema.js';
import * as graphicsShaderCode from '../shader-code.js';
import type { ShaderCode } from '../shader-code.js';

/** test helpers shared by the shader tests, not a test file itself */

export const ComputeShaders: [name: string, shaderCode: ShaderCode, schema: RenderableSchema][] = [
  ['background', background.SkyboxShaderCode, background.BackgroundSchema],
  ['bloom-luminosity', bloom.LuminosityShaderCode, bloom.LuminositySchema],
  ['bloom-blur', bloom.BlurShaderCode, bloom.BlurSchema],
  ['bloom-composite', bloom.CompositeShaderCode, bloom.CompositeSchema],
  ['cas', cas.CasShaderCode, cas.CasSchema],
  ['dof', dof.DofShaderCode, dof.DofSchema],
  ['dpoit-blend-back', dpoit.BlendBackDpoitShaderCode, dpoit.BlendBackDpoitSchema],
  ['dpoit-evaluate', dpoit.EvaluateDpoitShaderCode, dpoit.EvaluateDpoitSchema],
  ['fxaa', fxaa.FxaaShaderCode, fxaa.FxaaSchema],
  ['illumination-compose', illumination.ComposeShaderCode, illumination.ComposeSchema],
  [
    'illumination-multi-sample-compose',
    illumination.MultiSampleComposeShaderCode,
    illumination.MultiSampleComposeSchema,
  ],
  ['marking-dim', marking.DimShaderCode, marking.DimSchema],
  ['marking-overlay', marking.OverlayShaderCode, marking.OverlaySchema],
  ['multi-sample-compose', multiSample.ComposeShaderCode, multiSample.ComposeSchema],
  ['outlines', outline.OutlinesShaderCode, outline.OutlinesSchema],
  ['postprocessing', postprocessing.PostprocessingShaderCode, postprocessing.PostprocessingSchema],
  ['shadows', shadow.ShadowsShaderCode, shadow.ShadowsSchema],
  ['smaa-edges', smaa.EdgesShaderCode, smaa.EdgesSchema],
  ['smaa-weights', smaa.WeightsShaderCode, smaa.WeightsSchema],
  ['smaa-blend', smaa.BlendShaderCode, smaa.BlendSchema],
  ['ssao', ssao.SsaoShaderCode, ssao.SsaoSchema],
  ['ssao-blur', ssao.SsaoBlurShaderCode, ssao.SsaoBlurSchema],
  ['tracing-trace', tracing.TraceShaderCode, tracing.TraceSchema],
  ['wboit-evaluate', wboit.EvaluateWboitShaderCode, wboit.EvaluateWboitSchema],
  ['volume-composite', volume.VolumeCompositeShaderCode, volume.VolumeCompositeSchema],
  ['gaussian-density', gaussianDensity.GaussianDensityShaderCode, gaussianDensity.GaussianDensitySchema],
  ['color-accumulate', colorSmoothing.ColorAccumulateShaderCode, colorSmoothing.ColorAccumulateSchema],
  ['color-normalize', colorSmoothing.ColorNormalizeShaderCode, colorSmoothing.ColorNormalizeSchema],
  [
    'histogram-pyramid-reduction',
    histogramReduction.HistopyramidReductionShaderCode,
    histogramReduction.HistopyramidReductionSchema,
  ],
  ['histogram-pyramid-sum', histogramSum.HistopyramidSumShaderCode, histogramSum.HistopyramidSumSchema],
  ['hi-z', hiZ.HiZShaderCode, hiZ.HiZSchema],
  ['active-voxels', activeVoxels.ActiveVoxelsShaderCode, activeVoxels.ActiveVoxelsSchema],
  ['isosurface', isosurface.IsosurfaceShaderCode, isosurface.IsosurfaceSchema],
  ['copy', computeUtil.CopyShaderCode, computeUtil.CopySchema],
];

/** modules defining shader codes, by path relative to `packages/graphics/src` */
export const ShaderModules: { [path: string]: object } = {
  'canvas3d/passes/background.ts': background,
  'canvas3d/passes/bloom.ts': bloom,
  'canvas3d/passes/cas.ts': cas,
  'canvas3d/passes/dof.ts': dof,
  'canvas3d/passes/dpoit.ts': dpoit,
  'canvas3d/passes/fxaa.ts': fxaa,
  'canvas3d/passes/illumination.ts': illumination,
  'canvas3d/passes/marking.ts': marking,
  'canvas3d/passes/multi-sample.ts': multiSample,
  'canvas3d/passes/outline.ts': outline,
  'canvas3d/passes/postprocessing.ts': postprocessing,
  'canvas3d/passes/shadow.ts': shadow,
  'canvas3d/passes/smaa.ts': smaa,
  'canvas3d/passes/ssao.ts': ssao,
  'canvas3d/passes/tracing.ts': tracing,
  'canvas3d/passes/wboit.ts': wboit,
  'canvas3d/passes/volume.ts': volume,
  'geo/gaussian-density/gpu.ts': gaussianDensity,
  'geo/geometry/texture-mesh/color-smoothing.ts': colorSmoothing,
  'gl/webgl/compute/histogram-pyramid/reduction.ts': histogramReduction,
  'gl/webgl/compute/histogram-pyramid/sum.ts': histogramSum,
  'gl/compute/hi-z.ts': hiZ,
  'gl/webgl/compute/marching-cubes/active-voxels.ts': activeVoxels,
  'gl/webgl/compute/marching-cubes/isosurface.ts': isosurface,
  'gl/compute/util.ts': computeUtil,
  'gl/shader-code.ts': graphicsShaderCode,
};

/** example values for number defines, others are 1 */
const NumberDefines: { [k: string]: number } = {
  dNSamples: 32,
  dLevels: 3,
  dOcclusionKernelSize: 15,
  dMaxSearchSteps: 16,
  dEdgeThreshold: 0.1,
  dLightCount: 2,
  dOutlineScale: 2,
  dSteps: 4,
};

export type DefineValueMap = { [k: string]: string | boolean | number };

/**
 * Base variant (booleans on, first string option) and one variant per other string option or
 * disabled boolean; the declarations only depend on single defines (`usedIf`, `dimension`)
 */
export function getComputeDefineVariants(schema: RenderableSchema): DefineValueMap[] {
  const base: DefineValueMap = { dRenderVariant: 'compute' };
  for (const [k, spec] of Object.entries(schema)) {
    if (spec.type !== 'define') continue;
    if (spec.kind === 'string') {
      if (!spec.options?.length) throw new Error(`string define '${k}' has no options`);
      base[k] = spec.options[0];
    } else if (spec.kind === 'boolean') base[k] = true;
    else base[k] = NumberDefines[k] ?? 1;
  }
  const variants = [base];
  for (const [k, spec] of Object.entries(schema)) {
    if (spec.type !== 'define') continue;
    if (spec.kind === 'string') {
      for (const option of spec.options!.slice(1)) variants.push({ ...base, [k]: option });
    } else if (spec.kind === 'boolean') {
      variants.push({ ...base, [k]: false });
    }
  }
  return variants;
}

/** prefix defines of a program on a minimum-spec context, see `addShaderDefines` */
export function getComputePrefixDefines(shaderCode: ShaderCode, isWebGL2: boolean): string[] {
  const { drawBuffers, fragDepth, shaderTextureLod } = shaderCode.extensions;
  if (!isWebGL2) return drawBuffers === 'required' ? ['emulateVertexId', 'requiredDrawBuffers'] : ['emulateVertexId'];
  const defines = ['depthTextureSupport'];
  if (drawBuffers) defines.push('requiredDrawBuffers');
  if (fragDepth) defines.push('enabledFragDepth');
  if (shaderTextureLod) defines.push('enabledShaderTextureLod');
  return defines;
}
