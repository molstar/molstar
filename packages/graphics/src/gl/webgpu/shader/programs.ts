/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import type { ProgramVariant } from '../../render-item.js';
import type { DefineValues } from '../../shader-code.js';
import { copy_wgsl } from './copy.wgsl.js';
import { cylinders_wgsl } from './cylinders.wgsl.js';
import { dim_wgsl, edge_wgsl, overlay_wgsl } from './marking.wgsl.js';
import { evaluate_wboit_wgsl } from './wboit.wgsl.js';
import { blend_back_dpoit_wgsl, evaluate_dpoit_wgsl } from './dpoit.wgsl.js';
import { bloom_blur_wgsl, bloom_composite_wgsl, bloom_luminosity_wgsl } from './bloom.wgsl.js';
import { lines_wgsl } from './lines.wgsl.js';
import { mesh_wgsl } from './mesh.wgsl.js';
import { multi_sample_compose_wgsl } from './multi-sample.wgsl.js';
import { points_wgsl } from './points.wgsl.js';
import { background_wgsl } from './background.wgsl.js';
import { cas_wgsl } from './cas.wgsl.js';
import { direct_volume_wgsl } from './direct-volume.wgsl.js';
import { dof_wgsl } from './dof.wgsl.js';
import { illumination_compose_wgsl, trace_wgsl } from './illumination.wgsl.js';
import { image_wgsl } from './image.wgsl.js';
import { fxaa_wgsl } from './fxaa.wgsl.js';
import { outlines_wgsl } from './outlines.wgsl.js';
import { postprocessing_wgsl } from './postprocessing.wgsl.js';
import { shadows_wgsl } from './shadows.wgsl.js';
import { ssao_blur_wgsl } from './ssao-blur.wgsl.js';
import { smaa_blend_wgsl, smaa_edges_wgsl, smaa_weights_wgsl } from './smaa.wgsl.js';
import { ssao_wgsl } from './ssao.wgsl.js';
import { spheres_wgsl } from './spheres.wgsl.js';
import { text_wgsl } from './text.wgsl.js';

/** A WGSL program, by the name of the `ShaderCode` it replaces */
export interface WgslProgram {
  readonly source: string;
  /** render variants the source implements */
  readonly variants: readonly ProgramVariant[];
  /** the first feature selected by `defines` that the source doesn't implement */
  readonly getUnsupported?: (defines: DefineValues) => string | undefined;
}

function value(defines: DefineValues, name: string) {
  return defines[name]?.ref.value as unknown;
}

/** Render variants of the graphics programs, `color` with the transparency it is drawn with */
const PortedVariants = ['colorBlended', 'colorWboit', 'colorDpoit', 'pick', 'depth', 'marking', 'emissive', 'tracing'];

/** Data types other than these per-instance ones need volume or segment support */
const InstanceTypes = new Set(['uniform', 'instance', 'group', 'groupInstance', 'vertex', 'vertexInstance']);

/** Features the ported color variants share and don't implement yet */
function getUnsupported(
  geometryTypes: string | string[],
  defines: DefineValues,
  options: { solidInterior?: boolean } = {},
) {
  const variant = value(defines, 'dRenderVariant');
  if (!PortedVariants.includes(variant as string)) {
    return `render variant '${variant}'`;
  }
  const geometryType = value(defines, 'dGeometryType') as string;
  if (!([] as string[]).concat(geometryTypes).includes(geometryType)) return `geometry type '${geometryType}'`;
  const color = value(defines, 'dColorType') as string;
  if (!InstanceTypes.has(color)) return `color type '${color}'`;
  if (value(defines, 'dUsePalette')) return 'palettes';
  for (const type of ['dOverpaintType', 'dTransparencyType', 'dEmissiveType', 'dSubstanceType', 'dMarkerType']) {
    const v = value(defines, type) as string | undefined;
    if (v !== undefined && !InstanceTypes.has(v)) return `${type} '${v}'`;
  }
  if (value(defines, 'dSolidInterior') && !options.solidInterior) return 'solid interiors';
  if (value(defines, 'dSegmented')) return 'merged render items';
  if (value(defines, 'dWiggle')) return 'wiggle textures';
}

/** Features of images not implemented yet */
function getImageUnsupported(defines: DefineValues) {
  const variant = value(defines, 'dRenderVariant') as string;
  if (!PortedVariants.includes(variant)) {
    return `render variant '${variant}'`;
  }
}

/** Features of direct volumes not implemented yet; only the color variants are drawn, as on WebGL */
function getDirectVolumeUnsupported(defines: DefineValues) {
  const variant = value(defines, 'dRenderVariant') as string;
  if (!['colorBlended', 'colorWboit', 'colorDpoit'].includes(variant)) return `render variant '${variant}'`;
  if (value(defines, 'dGridTexType') !== '3d') return `grid texture type '${value(defines, 'dGridTexType')}'`;
  const color = value(defines, 'dColorType') as string;
  if (!InstanceTypes.has(color) && color !== 'direct') return `color type '${color}'`;
}

const Programs: { [name: string]: WgslProgram } = {
  mesh: {
    source: mesh_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: (d) => getUnsupported(['mesh', 'textureMesh'], d, { solidInterior: true }),
  },
  spheres: {
    source: spheres_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: (d) => getUnsupported('spheres', d, { solidInterior: true }),
  },
  cylinders: {
    source: cylinders_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: (d) => getUnsupported('cylinders', d, { solidInterior: true }),
  },
  copy: {
    source: copy_wgsl,
    variants: ['compute'],
    // packing binds depth textures as color, only needed without depth textures, which WebGPU always has
    getUnsupported: (d) => (value(d, 'dPackDepth') ? 'packed depth copies' : undefined),
  },
  ssao: { source: ssao_wgsl, variants: ['compute'] },
  edge: { source: edge_wgsl, variants: ['compute'] },
  'evaluate-wboit': { source: evaluate_wboit_wgsl, variants: ['compute'] },
  'Bloom Luminosity': { source: bloom_luminosity_wgsl, variants: ['compute'] },
  'Bloom Blur': { source: bloom_blur_wgsl, variants: ['compute'] },
  'Bloom Composite': { source: bloom_composite_wgsl, variants: ['compute'] },
  'evaluate-dpoit': { source: evaluate_dpoit_wgsl, variants: ['compute'] },
  'blend-back-dpoit': { source: blend_back_dpoit_wgsl, variants: ['compute'] },
  'multi-sample-compose': { source: multi_sample_compose_wgsl, variants: ['compute'] },
  trace: { source: trace_wgsl, variants: ['compute'] },
  // a weighted texel, as the multi-sample compose
  'trace-accumulate': { source: multi_sample_compose_wgsl, variants: ['compute'] },
  'illumination-compose': { source: illumination_compose_wgsl, variants: ['compute'] },
  'illumination-multi-sample-compose': { source: multi_sample_compose_wgsl, variants: ['compute'] },
  'smaa-edges': { source: smaa_edges_wgsl, variants: ['compute'] },
  'smaa-weights': { source: smaa_weights_wgsl, variants: ['compute'] },
  'smaa-blend': { source: smaa_blend_wgsl, variants: ['compute'] },
  dim: { source: dim_wgsl, variants: ['compute'] },
  overlay: { source: overlay_wgsl, variants: ['compute'] },
  ssao_blur: { source: ssao_blur_wgsl, variants: ['compute'] },
  postprocessing: { source: postprocessing_wgsl, variants: ['compute'] },
  outlines: { source: outlines_wgsl, variants: ['compute'] },
  shadows: { source: shadows_wgsl, variants: ['compute'] },
  dof: { source: dof_wgsl, variants: ['compute'] },
  fxaa: { source: fxaa_wgsl, variants: ['compute'] },
  cas: { source: cas_wgsl, variants: ['compute'] },
  background: { source: background_wgsl, variants: ['compute'] },
  lines: {
    source: lines_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: (d) => getUnsupported('lines', d),
  },
  image: {
    source: image_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: getImageUnsupported,
  },
  'direct-volume': { source: direct_volume_wgsl, variants: ['color'], getUnsupported: getDirectVolumeUnsupported },
  text: {
    source: text_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: (d) => getUnsupported('text', d),
  },
  points: {
    source: points_wgsl,
    variants: ['color', 'pick', 'depth', 'marking', 'emissive', 'tracing'],
    getUnsupported: (d) => getUnsupported('points', d),
  },
};

export function getWgslProgram(name: string): WgslProgram | undefined {
  return Programs[name];
}
