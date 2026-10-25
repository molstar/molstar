/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { getAttributeLayout } from '../attribute-layout.js';
import { MeshSchema } from '../renderable/mesh.js';
import {
  getDefineLookup,
  GlobalDefineSchema,
  GlobalTextureSchema,
  GlobalUniformSchema,
  InternalSchema,
} from '../renderable/schema.js';
import type { DefineValues } from '../shader-code.js';
import { preprocessShader } from '../shader-preprocessor.js';
import { getBindingLayout } from '../uniform-layout.js';
import {
  getOverrideConstants,
  getTextureSampleType,
  getWgslDefineLookup,
  getWgslDefines,
  isFilteredTexture,
} from '../webgpu/program.js';
import { copy_wgsl } from '../webgpu/shader/copy.wgsl.js';
import { postprocessing_wgsl } from '../webgpu/shader/postprocessing.wgsl.js';
import { ssao_blur_wgsl } from '../webgpu/shader/ssao-blur.wgsl.js';
import { ssao_wgsl } from '../webgpu/shader/ssao.wgsl.js';
import { evaluate_wboit_wgsl } from '../webgpu/shader/wboit.wgsl.js';
import { blend_back_dpoit_wgsl, evaluate_dpoit_wgsl } from '../webgpu/shader/dpoit.wgsl.js';
import { lines_wgsl } from '../webgpu/shader/lines.wgsl.js';
import { text_wgsl } from '../webgpu/shader/text.wgsl.js';
import { multi_sample_compose_wgsl } from '../webgpu/shader/multi-sample.wgsl.js';
import { illumination_compose_wgsl, trace_wgsl } from '../webgpu/shader/illumination.wgsl.js';
import { smaa_blend_wgsl, smaa_edges_wgsl, smaa_weights_wgsl } from '../webgpu/shader/smaa.wgsl.js';
import { dim_wgsl, overlay_wgsl } from '../webgpu/shader/marking.wgsl.js';
import { bloom_blur_wgsl, bloom_composite_wgsl, bloom_luminosity_wgsl } from '../webgpu/shader/bloom.wgsl.js';
import { outlines_wgsl } from '../webgpu/shader/outlines.wgsl.js';
import { shadows_wgsl } from '../webgpu/shader/shadows.wgsl.js';
import { background_wgsl } from '../webgpu/shader/background.wgsl.js';
import { cas_wgsl } from '../webgpu/shader/cas.wgsl.js';
import { dof_wgsl } from '../webgpu/shader/dof.wgsl.js';
import { fxaa_wgsl } from '../webgpu/shader/fxaa.wgsl.js';
import { direct_volume_wgsl } from '../webgpu/shader/direct-volume.wgsl.js';
import { image_wgsl } from '../webgpu/shader/image.wgsl.js';
import { mesh_wgsl } from '../webgpu/shader/mesh.wgsl.js';
import { CopySchema } from '../compute/util.js';
import { points_wgsl } from '../webgpu/shader/points.wgsl.js';
import { spheres_wgsl } from '../webgpu/shader/spheres.wgsl.js';
import { cylinders_wgsl } from '../webgpu/shader/cylinders.wgsl.js';
import { getWgslProgram, resolveWgslSource, type WgslProgram } from '../webgpu/shader/programs.js';
import {
  getWgslDeclarations,
  getWgslSamplerBinding,
  getWgslTextureDeclarations,
  getWgslTextureType,
  getWgslStorageDeclarations,
  getWgslVertexInputDeclaration,
} from '../webgpu/wgsl-declarations.js';
import { getWgslDrawBinding, getWgslStorageBindings } from '../webgpu/uniform-declarations.js';
import { StorageBufferSpec } from '../compute.js';

const MeshFullSchema = {
  ...GlobalUniformSchema,
  ...GlobalTextureSchema,
  ...GlobalDefineSchema,
  ...InternalSchema,
  ...MeshSchema,
};

function defineValues(values: { [k: string]: string | number | boolean }): DefineValues {
  const out: DefineValues = {};
  for (const k of Object.keys(values)) out[k] = ValueCell.create(values[k] as any);
  return out;
}

describe('wgsl declarations', () => {
  it('texture types', () => {
    expect(getWgslTextureType({ dimension: '2d', sampleType: 'float' })).toBe('texture_2d<f32>');
    expect(getWgslTextureType({ dimension: '2d', sampleType: 'unfilterable-float' })).toBe('texture_2d<f32>');
    expect(getWgslTextureType({ dimension: '3d', sampleType: 'float' })).toBe('texture_3d<f32>');
    expect(getWgslTextureType({ dimension: 'cube', sampleType: 'float' })).toBe('texture_cube<f32>');
    expect(getWgslTextureType({ dimension: '2d', sampleType: 'sint' })).toBe('texture_2d<i32>');
    expect(getWgslTextureType({ dimension: '2d', sampleType: 'uint' })).toBe('texture_2d<u32>');
    expect(getWgslTextureType({ dimension: '2d', sampleType: 'depth' })).toBe('texture_depth_2d');
    expect(() => getWgslTextureType({ dimension: '3d', sampleType: 'depth' })).toThrow();
  });

  it('textures and samplers use the binding layout', () => {
    const layout = getBindingLayout(MeshFullSchema, getDefineLookup([]));
    const code = getWgslTextureDeclarations(layout);
    for (const t of layout.object.textures) {
      expect(code).toContain(`@group(2) @binding(${t.binding}) var ${t.name}: `);
    }
    for (const t of layout.global.textures) {
      expect(code).toContain(`@group(0) @binding(${t.binding}) var ${t.name}: `);
    }
    const samplerBinding = getWgslSamplerBinding(layout);
    expect(samplerBinding).toBe(layout.global.textures.length + 1);
    expect(code).toContain(`@group(0) @binding(${samplerBinding}) var sNearest: sampler;`);
    expect(code).toContain(`@group(0) @binding(${samplerBinding + 2}) var sLinearMipmap: sampler;`);

    // the theme layers (`DataSpec`) are no textures on WebGPU: read from the data pages
    expect(layout.object.textures.find((t) => t.name === 'tColor')).toBeUndefined();
    const gridLayout = getBindingLayout(MeshFullSchema, getDefineLookup(['dColorType_grid']));
    const tColorGrid = gridLayout.object.textures.find((t) => t.name === 'tColorGrid')!;
    const referenced = getWgslTextureDeclarations(gridLayout, (n) => n === 'tColorGrid' || n === 'sLinear');
    expect(referenced.split('\n')).toEqual([
      `@group(2) @binding(${tColorGrid.binding}) var tColorGrid: texture_2d<f32>;`,
      `@group(0) @binding(${samplerBinding + 1}) var sLinear: sampler;`,
    ]);
  });

  it('vertex inputs use the attribute layout', () => {
    const lookup = getWgslDefineLookup(new Map());
    const code = getWgslVertexInputDeclaration(MeshFullSchema, lookup);
    const { byName } = getAttributeLayout(MeshFullSchema);
    expect(code).toContain('@builtin(vertex_index) vertexIndex: u32,');
    expect(code).toContain('@builtin(instance_index) instanceIndex: u32,');
    expect(code).toContain(`@location(${byName.aPosition.location}) aPosition: vec3<f32>,`);
    expect(code).toContain(`@location(${byName.aGroup.location}) aGroup: f32,`);
    for (let i = 0; i < 4; ++i) {
      expect(code).toContain(`@location(${byName.aTransform.location + i}) aTransform${i}: vec4<f32>,`);
    }
    // vertex ids are builtins
    expect(code).not.toContain('aVertex');
  });
});

describe('wgsl defines', () => {
  it('flattens define values', () => {
    const defines = getWgslDefines(defineValues({ dA: 'x', dB: true, dC: false, dLightCount: 3 }));
    expect([...defines]).toEqual([
      ['dA_x', ''],
      ['dB', ''],
      ['dLightCount', '3'],
    ]);
  });

  it('override constants', () => {
    const source = 'override oLightCount: i32 = 1;\noverride oOther: f32;';
    expect(getOverrideConstants(source, defineValues({ dLightCount: 2 }))).toEqual({ oLightCount: 2 });
  });

  it('override constants per entry point', () => {
    const source = [
      'override oLightCount: i32 = 1;',
      'override oSteps: i32 = 4;',
      'override oClipObjectCount: i32 = 1;',
      'const kSteps = oSteps * 2; // oClipObjectCount',
      'struct Light { count: i32 }',
      'fn light() -> i32 { if (true) { return oLightCount; } return 0; }',
      'fn steps() -> i32 { return kSteps; }',
      '@vertex fn vertexMain() -> @builtin(position) vec4f { let n = steps(); return vec4f(f32(n)); }',
      '@fragment fn fragmentMain() -> @location(0) vec4f { let l = Light(light()); return vec4f(f32(l.count)); }',
    ].join('\n');
    const values = defineValues({ dLightCount: 2, dSteps: 8, dClipObjectCount: 3 });
    expect(getOverrideConstants(source, values, 'vertexMain')).toEqual({ oSteps: 8 });
    expect(getOverrideConstants(source, values, 'fragmentMain')).toEqual({ oLightCount: 2 });
    expect(getOverrideConstants(source, values)).toEqual({ oLightCount: 2, oSteps: 8, oClipObjectCount: 3 });
  });
});

describe('mesh wgsl', () => {
  const base = {
    dRenderVariant: 'colorBlended',
    dGeometryType: 'mesh',
    dColorType: 'texture',
    dLightCount: 1,
    dClipObjectCount: 0,
    dColorMarker: true,
  };
  const cases: { [k: string]: string | number | boolean }[] = [
    base,
    { ...base, dColorType: 'uniform', dFlatShaded: true, dIgnoreLight: true },
    { ...base, dColorType: 'texture', dCelShaded: true, dXrayShaded: 'inverted' },
    { ...base, dColorType: 'grid' },
    {
      ...base,
      dColorType: 'grid',
      dOverpaint: true,
      dOverpaintType: 'volumeInstance',
      dTransparency: true,
      dTransparencyType: 'volumeInstance',
      dEmissive: true,
      dEmissiveType: 'volumeInstance',
      dSubstance: true,
      dSubstanceType: 'volumeInstance',
    },
    {
      ...base,
      dOverpaint: true,
      dOverpaintType: 'texture',
      dTransparency: true,
      dTransparencyType: 'texture',
      dEmissive: true,
      dEmissiveType: 'vertexInstance',
      dSubstance: true,
      dSubstanceType: 'texture',
      dTransparentBackfaces: 'opaque',
      dFlipSided: true,
    },
    { ...base, dColorMarker: false, dTransparentBackfaces: 'off', dXrayShaded: 'on' },
    { ...base, dRenderVariant: 'pick', dTransparency: true, dTransparencyType: 'texture', dXrayShaded: 'on' },
    {
      ...base,
      dRenderVariant: 'marking',
      dTransparentBackfaces: 'off',
      dTransparency: true,
      dTransparencyType: 'texture',
    },
    { ...base, dRenderVariant: 'colorWboit', dTransparentBackfaces: 'off' },
    { ...base, dRenderVariant: 'colorDpoit', dTransparentBackfaces: 'off', dXrayShaded: 'on' },
    {
      ...base,
      dRenderVariant: 'emissive',
      dEmissive: true,
      dEmissiveType: 'texture',
      dOverpaint: true,
      dOverpaintType: 'texture',
      dTransparency: true,
      dTransparencyType: 'texture',
      dXrayShaded: 'on',
    },
    { ...base, dRenderVariant: 'emissive', dColorType: 'uniform' },
    { ...base, dRenderVariant: 'tracing', dTransparentBackfaces: 'opaque', dEmissive: true, dEmissiveType: 'texture' },
    { ...base, dRenderVariant: 'tracing', dSolidInterior: true, dIgnoreLight: true },
    { ...base, dGeometryType: 'textureMesh', dVaryingGroup: true },
    { ...base, dGeometryType: 'textureMesh', dRenderVariant: 'pick', dFlatShaded: true },
    { ...base, dRenderVariant: 'depth', dTransparency: true, dTransparencyType: 'texture', dVaryingGroup: true },
    { ...base, dRenderVariant: 'depth', dTransparentBackfaces: 'off', dXrayShaded: 'on' },
    { ...base, dSolidInterior: true, dFlatShaded: true },
    {
      ...base,
      dRenderVariant: 'pick',
      dSolidInterior: true,
      dClipObjects: true,
      dClipObjectCount: 2,
      dClipVariant: 'pixel',
      dClipping: true,
    },
    {
      ...base,
      dRenderVariant: 'depth',
      dSolidInterior: true,
      dClipObjects: true,
      dClipObjectCount: 1,
      dClipVariant: 'pixel',
    },
    { ...base, dClipObjects: true, dClipObjectCount: 2, dClipVariant: 'pixel', dClipping: true },
    { ...base, dClipObjects: true, dClipObjectCount: 1, dClipVariant: 'instance', dClipping: true },
    {
      ...base,
      dRenderVariant: 'pick',
      dClipObjects: true,
      dClipObjectCount: 3,
      dClipVariant: 'pixel',
      dVaryingGroup: true,
      dClipping: true,
    },
  ];

  for (const values of cases) {
    it(`preprocesses ${JSON.stringify(values)}`, () => {
      const dv = defineValues(values);
      expect(getWgslProgram('mesh')!.getUnsupported!(dv)).toBeUndefined();
      const code = preprocessShader(mesh_wgsl, getWgslDefines(dv), { expandMacros: true });
      expect(code).toContain('fn vertexMain(');
      expect(code).toContain('fn fragmentMain(');
      expect(code).not.toMatch(/^\s*#/m);
    });
  }

  it('reports unsupported features', () => {
    const { getUnsupported } = getWgslProgram('mesh')!;
    expect(getUnsupported!(defineValues({ ...base, dUsePalette: true }))).toBeUndefined();
    expect(getUnsupported!(defineValues({ ...base, dWiggle: true }))).toBeUndefined();
    // grid data is read by the mesh programs only
    const spheres = getWgslProgram('spheres')!;
    expect(spheres.getUnsupported!(defineValues({ ...base, dGeometryType: 'spheres', dColorType: 'grid' }))).toBe(
      "color type 'grid'",
    );
  });

  it('reports unsupported direct volume features', () => {
    const { getUnsupported } = getWgslProgram('direct-volume')!;
    const dv = { ...base, dGeometryType: 'directVolume', dGridTexType: '3d' };
    expect(getUnsupported!(defineValues({ ...dv, dGridTexType: '2d' }))).toBe("grid texture type '2d'");
    expect(getUnsupported!(defineValues({ ...dv, dRenderVariant: 'pick' }))).toBeUndefined();
    expect(getUnsupported!(defineValues({ ...dv, dRenderVariant: 'tracing' }))).toBe("render variant 'tracing'");
  });
});

describe('quad wgsl', () => {
  const clip = { dClipObjects: true, dClipObjectCount: 2, dClipVariant: 'pixel', dClipping: true };
  const base = {
    dRenderVariant: 'colorBlended',
    dColorType: 'texture',
    dSizeType: 'uniform',
    dLightCount: 1,
    dClipObjectCount: 0,
    dColorMarker: true,
  };
  const cases: [string, WgslProgram['source'], { [k: string]: string | number | boolean }][] = [
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dSolidInterior: true }],
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dApproximate: true, dSizeType: 'texture' }],
    [
      'spheres',
      spheres_wgsl,
      { ...base, dGeometryType: 'spheres', dColorType: 'texture', dSizeType: 'texture', dXrayShaded: 'on' },
    ],
    ['points', points_wgsl, { ...base, dGeometryType: 'points', dPointStyle: 'circle', dSizeType: 'texture' }],
    [
      'points',
      points_wgsl,
      { ...base, dGeometryType: 'points', dPointStyle: 'fuzzy', dPointSizeAttenuation: true, dSizeType: 'texture' },
    ],
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dRenderVariant: 'pick' }],
    ['points', points_wgsl, { ...base, dGeometryType: 'points', dRenderVariant: 'pick', dPointStyle: 'circle' }],
    ['cylinders', cylinders_wgsl, { ...base, dGeometryType: 'cylinders', dDualColor: true }],
    [
      'cylinders',
      cylinders_wgsl,
      { ...base, dGeometryType: 'cylinders', dColorType: 'texture', dDualColor: true, dSizeType: 'texture' },
    ],
    ['cylinders', cylinders_wgsl, { ...base, dGeometryType: 'cylinders', dRenderVariant: 'pick', dXrayShaded: 'on' }],
    [
      'cylinders',
      cylinders_wgsl,
      { ...base, dGeometryType: 'cylinders', dRenderVariant: 'marking', dSolidInterior: true },
    ],
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dRenderVariant: 'marking' }],
    ['points', points_wgsl, { ...base, dGeometryType: 'points', dRenderVariant: 'marking', dPointStyle: 'fuzzy' }],
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dRenderVariant: 'colorWboit' }],
    ['cylinders', cylinders_wgsl, { ...base, dGeometryType: 'cylinders', dRenderVariant: 'colorDpoit' }],
    ['lines', lines_wgsl, { ...base, dGeometryType: 'lines', dLineSizeAttenuation: true }],
    ['lines', lines_wgsl, { ...base, dGeometryType: 'lines', dRenderVariant: 'pick' }],
    ['lines', lines_wgsl, { ...base, dGeometryType: 'lines', dRenderVariant: 'colorWboit' }],
    ['text', text_wgsl, { ...base, dGeometryType: 'text' }],
    ['text', text_wgsl, { ...base, dGeometryType: 'text', dRenderVariant: 'marking' }],
    ['text', text_wgsl, { ...base, dGeometryType: 'text', dRenderVariant: 'colorDpoit' }],
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dRenderVariant: 'emissive', dXrayShaded: 'on' }],
    [
      'cylinders',
      cylinders_wgsl,
      { ...base, dGeometryType: 'cylinders', dRenderVariant: 'emissive', dDualColor: true },
    ],
    ['points', points_wgsl, { ...base, dGeometryType: 'points', dRenderVariant: 'emissive', dPointStyle: 'circle' }],
    ['lines', lines_wgsl, { ...base, dGeometryType: 'lines', dRenderVariant: 'emissive' }],
    ['text', text_wgsl, { ...base, dGeometryType: 'text', dRenderVariant: 'emissive' }],
    ['spheres', spheres_wgsl, { ...base, dGeometryType: 'spheres', dRenderVariant: 'tracing' }],
    ['cylinders', cylinders_wgsl, { ...base, dGeometryType: 'cylinders', dRenderVariant: 'tracing' }],
    ['points', points_wgsl, { ...base, dGeometryType: 'points', dRenderVariant: 'tracing', dPointStyle: 'circle' }],
    ['lines', lines_wgsl, { ...base, dGeometryType: 'lines', dRenderVariant: 'tracing' }],
    ['text', text_wgsl, { ...base, dGeometryType: 'text', dRenderVariant: 'tracing' }],
    ['image', image_wgsl, { ...base, dGeometryType: 'image', dInterpolation: 'nearest', dRenderVariant: 'tracing' }],
    [
      'spheres',
      spheres_wgsl,
      { ...base, ...clip, dGeometryType: 'spheres', dRenderVariant: 'depth', dSolidInterior: true },
    ],
    ['cylinders', cylinders_wgsl, { ...base, dGeometryType: 'cylinders', dRenderVariant: 'depth', dXrayShaded: 'on' }],
    ['points', points_wgsl, { ...base, dGeometryType: 'points', dRenderVariant: 'depth', dPointStyle: 'fuzzy' }],
    ['lines', lines_wgsl, { ...base, dGeometryType: 'lines', dRenderVariant: 'depth' }],
    ['text', text_wgsl, { ...base, dGeometryType: 'text', dRenderVariant: 'depth' }],
    ['spheres', spheres_wgsl, { ...base, ...clip, dGeometryType: 'spheres', dSolidInterior: true }],
    ['spheres', spheres_wgsl, { ...base, ...clip, dGeometryType: 'spheres', dClipPrimitive: true }],
    ['spheres', spheres_wgsl, { ...base, ...clip, dGeometryType: 'spheres', dClipVariant: 'instance' }],
    ['spheres', spheres_wgsl, { ...base, ...clip, dGeometryType: 'spheres', dApproximate: true }],
    ['cylinders', cylinders_wgsl, { ...base, ...clip, dGeometryType: 'cylinders', dSolidInterior: true }],
    ['cylinders', cylinders_wgsl, { ...base, ...clip, dGeometryType: 'cylinders', dRenderVariant: 'pick' }],
    ['points', points_wgsl, { ...base, ...clip, dGeometryType: 'points', dPointStyle: 'square' }],
    ['lines', lines_wgsl, { ...base, ...clip, dGeometryType: 'lines' }],
    ['text', text_wgsl, { ...base, ...clip, dGeometryType: 'text', dClipVariant: 'instance' }],
    ['image', image_wgsl, { ...base, ...clip, dGeometryType: 'image', dInterpolation: 'nearest' }],
    [
      'image',
      image_wgsl,
      {
        ...base,
        dGeometryType: 'image',
        dInterpolation: 'mitchell',
        dUsePalette: true,
        dOverpaint: true,
        dOverpaintType: 'texture',
      },
    ],
    [
      'image',
      image_wgsl,
      {
        ...base,
        dGeometryType: 'image',
        dInterpolation: 'bspline',
        dRenderVariant: 'pick',
        dTransparency: true,
        dTransparencyType: 'texture',
      },
    ],
    ['image', image_wgsl, { ...base, dGeometryType: 'image', dInterpolation: 'catmulrom', dRenderVariant: 'marking' }],
    ['image', image_wgsl, { ...base, dGeometryType: 'image', dInterpolation: 'nearest', dRenderVariant: 'depth' }],
    [
      'direct-volume',
      direct_volume_wgsl,
      {
        ...base,
        ...clip,
        dGeometryType: 'directVolume',
        dGridTexType: '3d',
        dColorType: 'direct',
        dUsePalette: true,
        dAxisOrder: '012',
      },
    ],
    [
      'direct-volume',
      direct_volume_wgsl,
      {
        ...base,
        dGeometryType: 'directVolume',
        dGridTexType: '3d',
        dColorType: 'texture',
        dPackedGroup: true,
        dRenderVariant: 'colorWboit',
        dOverpaint: true,
        dOverpaintType: 'vertexInstance',
      },
    ],
    [
      'direct-volume',
      direct_volume_wgsl,
      {
        ...base,
        dGeometryType: 'directVolume',
        dGridTexType: '3d',
        dIgnoreLight: true,
        dAxisOrder: '210',
        dRenderVariant: 'colorDpoit',
      },
    ],
  ];

  for (const [name, source, values] of cases) {
    it(`preprocesses ${name} ${JSON.stringify(values)}`, () => {
      const dv = defineValues(values);
      expect(getWgslProgram(name)!.getUnsupported!(dv)).toBeUndefined();
      const code = preprocessShader(resolveWgslSource(source, dv), getWgslDefines(dv), { expandMacros: true });
      expect(code).toContain('fn vertexMain(');
      expect(code).toContain('fn fragmentMain(');
      expect(code).not.toMatch(/^\s*#/m);
    });
  }
});

describe('compute wgsl', () => {
  const cases: [string, WgslProgram['source'], { [k: string]: string | number | boolean }][] = [
    ['ssao', ssao_wgsl, { dIncludeTransparent: true, dNSamples: 32, dLevels: 3 }],
    [
      'ssao',
      ssao_wgsl,
      { dMultiScale: true, dPackedDepthLevels: true, dIllumination: true, dNSamples: 16, dLevels: 3 },
    ],
    ['ssao_blur', ssao_blur_wgsl, { dOrthographic: 1, dOcclusionKernelSize: 15 }],
    [
      'postprocessing',
      postprocessing_wgsl,
      { dOrthographic: 0, dOcclusionEnable: true, dBlendTransparency: true, dOcclusionIncludeTransparency: true },
    ],
    ['dim', dim_wgsl, {}],
    ['overlay', overlay_wgsl, { dMarkingShading: 'ssao' }],
    ['overlay', overlay_wgsl, { dMarkingShading: 'off' }],
    ['overlay', overlay_wgsl, { dMarkingShading: 'traced' }],
    [
      'trace',
      trace_wgsl,
      { dOrthographic: 0, dLightCount: 2, dShadowEnable: true, dThicknessMode: 'auto', dSteps: 32, dBounces: 4 },
    ],
    ['trace', trace_wgsl, { dOrthographic: 1, dLightCount: 0, dThicknessMode: 'fixed', dRefineSteps: 0 }],
    [
      'illumination-compose',
      illumination_compose_wgsl,
      { dOrthographic: 0, dDenoise: true, dBlendTransparency: true, dOcclusionEnable: true, dBloomEnable: true },
    ],
    [
      'illumination-compose',
      illumination_compose_wgsl,
      { dOrthographic: 1, dOutlineEnable: true, dOutlineScale: 2, dTransparentOutline: true, dBlendTransparency: true },
    ],
    ['illumination-multi-sample-compose', multi_sample_compose_wgsl, {}],
    ['multi-sample-compose', multi_sample_compose_wgsl, {}],
    ['smaa-edges', smaa_edges_wgsl, { dEdgeThreshold: 0.1 }],
    ['smaa-weights', smaa_weights_wgsl, { dMaxSearchSteps: 16 }],
    ['smaa-blend', smaa_blend_wgsl, {}],
    ['evaluate-wboit', evaluate_wboit_wgsl, {}],
    ['evaluate-dpoit', evaluate_dpoit_wgsl, {}],
    ['blend-back-dpoit', blend_back_dpoit_wgsl, {}],
    ['postprocessing', postprocessing_wgsl, { dOrthographic: 0, dBloomEnable: true }],
    ['Bloom Luminosity', bloom_luminosity_wgsl, { dMode: 'luminosity', dBlendTransparency: true }],
    ['Bloom Luminosity', bloom_luminosity_wgsl, { dMode: 'luminosity' }],
    ['Bloom Luminosity', bloom_luminosity_wgsl, { dMode: 'emissive', dBlendTransparency: true }],
    ['Bloom Blur', bloom_blur_wgsl, { dKernelRadius: 5 }],
    ['Bloom Composite', bloom_composite_wgsl, {}],
    [
      'postprocessing',
      postprocessing_wgsl,
      { dOrthographic: 1, dOutlineEnable: true, dOutlineScale: 2, dShadowEnable: true, dBlendTransparency: true },
    ],
    ['outlines', outlines_wgsl, { dOrthographic: 0, dTransparentOutline: true }],
    ['outlines', outlines_wgsl, { dOrthographic: 1 }],
    ['shadows', shadows_wgsl, { dOrthographic: 0, dSteps: 8, dLightCount: 2 }],
    ['shadows', shadows_wgsl, { dOrthographic: 1, dSteps: 8, dLightCount: 0 }],
    ['dof', dof_wgsl, { dOrthographic: 0, dBlurSize: 9 }],
    [
      'fxaa',
      fxaa_wgsl,
      { dEdgeThresholdMin: 0.0312, dEdgeThresholdMax: 0.063, dIterations: 12, dSubpixelQuality: 0.3 },
    ],
    ['cas', cas_wgsl, { dDenoise: true }],
    ['background', background_wgsl, { dVariant: 'skybox' }],
    ['background', background_wgsl, { dVariant: 'image' }],
    ['background', background_wgsl, { dVariant: 'horizontalGradient' }],
    ['background', background_wgsl, { dVariant: 'radialGradient' }],
  ];
  for (const [name, source, values] of cases) {
    it(`preprocesses ${name} ${JSON.stringify(values)}`, () => {
      const dv = defineValues({ dRenderVariant: 'compute', ...values });
      expect(getWgslProgram(name)!.getUnsupported?.(dv)).toBeUndefined();
      const code = preprocessShader(resolveWgslSource(source, dv), getWgslDefines(dv), { expandMacros: true });
      expect(code).toContain('fn fragmentMain(');
      expect(code).not.toMatch(/^\s*#/m);
    });
  }

  it('copy', () => {
    const dv = defineValues({ dRenderVariant: 'compute', dPackDepth: false });
    const program = getWgslProgram('copy')!;
    expect(program.variants).toEqual(['compute']);
    expect(program.getUnsupported!(dv)).toBeUndefined();
    expect(program.getUnsupported!(defineValues({ dRenderVariant: 'compute', dPackDepth: true }))).toBe(
      'packed depth copies',
    );

    const defines = getWgslDefines(dv);
    const code = preprocessShader(copy_wgsl, defines, { expandMacros: true });
    expect(code).not.toMatch(/^\s*#/m);
    const lookup = getWgslDefineLookup(defines);
    const declarations = getWgslDeclarations(CopySchema, getBindingLayout(CopySchema, lookup), lookup);
    expect(declarations.unsupported).toEqual([]);
    expect(declarations.code).toContain('var tColor: texture_2d<f32>;');
    expect(declarations.code).toMatch(/@location\(\d+\) aPosition: vec2<f32>,/);
    expect(declarations.code).toContain('uQuadScale: vec2<f32>,');
  });
});

describe('wgsl texture sample types', () => {
  it('filtered textures', () => {
    expect(isFilteredTexture('let c = textureSampleLevel(tColor, sLinear, uv, 0.0);', 'tColor')).toBe(true);
    expect(isFilteredTexture('let c = textureSample( tColor, sLinear, uv);', 'tColor')).toBe(true);
    expect(isFilteredTexture('let c = textureLoad(tColor, coords, 0);', 'tColor')).toBe(false);
    expect(isFilteredTexture('let c = textureSampleLevel(tColorGrid, sLinear, uv, 0.0);', 'tColor')).toBe(false);
  });

  it('textures read with textureLoad accept unfilterable and depth formats', () => {
    expect(getTextureSampleType('float', false)).toBe('unfilterable-float');
    expect(getTextureSampleType('float', true)).toBe('float');
    expect(getTextureSampleType('depth', false)).toBe('depth');
    expect(getTextureSampleType('sint', false)).toBe('sint');
  });
});

describe('wgsl storage buffers of render schemas', () => {
  const schema = {
    ...MeshFullSchema,
    sA: StorageBufferSpec('vec4<f32>', 'read', ['dUseA']),
    sB: StorageBufferSpec('u32'),
  };

  it('numbers them after the draw uniforms, the same with and without the define', () => {
    const draw = getWgslDrawBinding(schema);
    const withA = getWgslStorageBindings(schema, (n) => n === 'dUseA');
    const withoutA = getWgslStorageBindings(schema, () => false);
    expect(withA).toEqual([
      { name: 'sA', binding: draw + 1, element: 'vec4<f32>' },
      { name: 'sB', binding: draw + 2, element: 'u32' },
    ]);
    expect(withoutA).toEqual([{ name: 'sB', binding: draw + 2, element: 'u32' }]);
  });

  it('declares referenced ones as read-only arrays in the object group', () => {
    const code = getWgslStorageDeclarations(
      schema,
      () => true,
      (n) => n === 'sB',
    );
    expect(code).toBe(`@group(2) @binding(${getWgslDrawBinding(schema) + 2}) var<storage, read> sB: array<u32>;`);
  });

  it('rejects writable storage buffers', () => {
    const writable = { ...MeshFullSchema, sC: StorageBufferSpec('u32', 'read-write') };
    expect(() => getWgslStorageBindings(writable, () => true)).toThrow(/read-only/);
  });
});
