/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { getAttributeLayout } from '../attribute-layout.js';
import {
  GlobalDefineSchema,
  GlobalTextureSchema,
  GlobalUniformSchema,
  getDefineLookup,
  InternalSchema,
  isUsed,
  type RenderableSchema,
  SegmentSchema,
} from '../renderable/schema.js';
import { getBindingLayout, getUniformBlockLayout, UniformBlocks } from '../uniform-layout.js';
import { ComputeShaders, getComputeDefineVariants, getComputePrefixDefines } from './compute-shaders.js';
import { GraphicsShaders, getProgramDefines, preprocessGlsl, preprocessProgram } from './graphics-shaders.js';
import { getGlslProgram } from '../webgl/programs.js';

/**
 * Everything must work with these; exceedances need a documented fallback and are listed
 * explicitly below.
 */
const WebGL2Minimum = {
  vertexAttribs: 16,
  uniformBlockSize: 16384,
  uniformBlocksPerStage: 12,
  texturesPerStage: 16,
  varyingVectors: 15,
};
const WebGPUDefault = {
  vertexAttributes: 16,
  vertexBuffers: 8,
  uniformBufferBindingSize: 65536,
  bindGroups: 4,
  sampledTexturesPerStage: 16,
  interStageVariables: 16,
};

/** one buffer per attribute; fallback (open): interleaving or storage buffers + vertex pulling */
const KnownVertexBufferExceedances = new Set(['cylinders', 'cylinders (merged)', 'lines (merged)', 'text (merged)']);

/** all textures without `usedIf` are bound; fallback (open): `usedIf` for the BaseSchema feature textures */
const KnownBindingLayoutExceedances = new Set<string>();

/**
 * WGSL has no automatic varying packing, every inter-stage variable takes a location (a `mat4` four);
 * fallback (open): pack scalar varyings into `vec4`s in the WGSL shaders
 */
const KnownInterStageExceedances = new Set(
  ['colorBlended', 'colorWboit', 'colorDpoit', 'tracing'].map((v) => `cylinders ${v}`),
);

function getSchema(schema: RenderableSchema, merged: boolean): RenderableSchema {
  return {
    ...GlobalUniformSchema,
    ...GlobalTextureSchema,
    ...GlobalDefineSchema,
    ...InternalSchema,
    ...schema,
    ...(merged ? SegmentSchema : {}),
  };
}

describe('minimum limits: attributes', () => {
  for (const [name, , schema] of GraphicsShaders) {
    for (const merged of [false, true]) {
      const id = merged ? `${name} (merged)` : name;
      const layout = getAttributeLayout(getSchema(schema, merged));
      it(`locations: ${id}`, () => {
        expect(layout.size).toBeLessThanOrEqual(Math.min(WebGL2Minimum.vertexAttribs, WebGPUDefault.vertexAttributes));
      });
      it(`webgpu vertex buffers: ${id}`, () => {
        // WebGPU has `vertex_index`, no `aVertex` emulation
        const count = layout.entries.filter((e) => e.name !== 'aVertex').length;
        expect(count > WebGPUDefault.vertexBuffers).toBe(KnownVertexBufferExceedances.has(id));
      });
    }
  }
});

describe('minimum limits: uniform blocks', () => {
  it('bind groups', () => {
    expect(UniformBlocks.length).toBeLessThanOrEqual(WebGPUDefault.bindGroups);
    expect(UniformBlocks.length).toBeLessThanOrEqual(WebGL2Minimum.uniformBlocksPerStage);
  });
  for (const [name, , schema] of GraphicsShaders) {
    it(`block sizes: ${name}`, () => {
      const s = getSchema(schema, true);
      for (const block of UniformBlocks) {
        const size = getUniformBlockLayout(block, s)?.size ?? 0;
        expect(size).toBeLessThanOrEqual(
          Math.min(WebGL2Minimum.uniformBlockSize, WebGPUDefault.uniformBufferBindingSize),
        );
      }
    });
  }
});

/** worst case: every feature that adds a texture or varying, filtered to the defines a schema has */
const FeatureDefines: { [k: string]: string | boolean | number } = {
  dLightCount: 1,
  dColorMarker: true,
  dClipObjects: true,
  dClipVariant: 'pixel',
  dColorType: 'texture',
  dUsePalette: true,
  dSizeType: 'texture',
  dOverpaint: true,
  dOverpaintType: 'texture',
  dTransparency: true,
  dTransparencyType: 'texture',
  dEmissive: true,
  dEmissiveType: 'texture',
  dSubstance: true,
  dSubstanceType: 'texture',
  dClipping: true,
  dWiggle: true,
  dSegmented: true,
  dXrayShaded: 'on',
  dTransparentBackfaces: 'on',
  dSolidInterior: true,
  dDualColor: true,
  dVaryingGroup: true,
  dGridTexType: '3d',
  dInterpolation: 'catmulrom',
  dPointStyle: 'fuzzy',
};

const RenderVariants = ['colorBlended', 'colorWboit', 'colorDpoit', 'pick', 'depth', 'marking', 'emissive', 'tracing'];

/** WebGL2 core features, no extensions */
const WebGL2PrefixDefines = [
  'requiredDrawBuffers',
  'enabledFragDepth',
  'enabledShaderTextureLod',
  'depthTextureSupport',
];

/**
 * schema textures and hand-declared samplers that are referenced (declarations are generated); data arrays read by
 * their accessors (`read_<name>`), a data texture each on WebGL
 */
function getActiveSamplers(source: string, schema: RenderableSchema) {
  const declaration = /uniform\s+(?:(?:lowp|mediump|highp)\s+)?\w*sampler\w+\s+(\w+)\s*;/g;
  const names = new Set([
    ...Object.keys(schema).filter((k) => schema[k].type === 'texture' || schema[k].type === 'data'),
    ...Array.from(source.matchAll(declaration), (m) => m[1]),
  ]);
  const code = source.replace(declaration, '');
  return [...names].filter((name) => new RegExp(`\\b(?:read_)?${name}\\b`).test(code));
}

function getVaryingTypes(vertSource: string) {
  return [
    ...vertSource.matchAll(/^\s*(?:flat\s+)?(?:varying|out)\s+(?:(?:lowp|mediump|highp)\s+)?(\w+)\s+\w+\s*;/gm),
  ].map((m) => m[1]);
}

/** inter-stage locations without packing (WGSL) */
function getInterStageLocations(types: string[]) {
  return types.reduce((count, type) => count + (type === 'mat4' ? 4 : type === 'mat3' ? 3 : 1), 0);
}

/** vec4 rows after GLSL ES 3.00 varying packing (simplified, appendix A.7) */
function getPackedVaryingRows(types: string[]) {
  let full = 0,
    threeColumns = 0,
    twoColumns = 0,
    singles = 0;
  for (const type of types) {
    if (type === 'mat4') full += 4;
    else if (type === 'mat3') threeColumns += 3;
    else if (/^[iu]?vec4$/.test(type)) full += 1;
    else if (/^[iu]?vec3$/.test(type)) threeColumns += 1;
    else if (/^[iu]?vec2$/.test(type)) twoColumns += 1;
    else singles += 1;
  }
  // singles fill the free column of three-column rows, then the free half of an odd two-column row
  singles = Math.max(0, singles - threeColumns);
  const twoColumnRows = Math.ceil(twoColumns / 2);
  if (twoColumns % 2 === 1) singles = Math.max(0, singles - 2);
  return full + threeColumns + twoColumnRows + Math.ceil(singles / 4);
}

/** `v*` identifiers used in a stage but not declared there (as varying or local) */
function getUndeclaredVaryings(source: string) {
  const code = source.replace(/\/\/.*$/gm, '');
  const declared = new Set(
    [...code.matchAll(/\b(?:(?:lowp|mediump|highp)\s+)?(?:float|int|[iu]?vec[234]|mat[234])\s+(v[A-Z]\w*)/g)].map(
      (m) => m[1],
    ),
  );
  return [...new Set(code.match(/\bv[A-Z]\w*/g) ?? [])].filter((n) => !declared.has(n));
}

describe('minimum limits: per program', () => {
  it('reference preprocessor', () => {
    const source = [
      '#define dA',
      '#if defined(dA) && !defined(dB)',
      'a',
      '#elif dN > 1',
      'b',
      '#else',
      'c',
      '#endif',
      '#ifdef dB',
      '#define dC',
      '#endif',
      '#ifndef dC',
      'd',
      '#endif',
    ].join('\n');
    expect(preprocessGlsl(source, new Map()).split('\n')).toEqual(['a', 'd']);
    expect(preprocessGlsl(source.replace('#define dA\n', ''), new Map([['dN', '2']])).split('\n')).toEqual(['b', 'd']);
    const expression = ['#if 1 + 2 * 3 == 7 && !(0 || dN < 2) && 10 - 4 - 3 == 3', 'e', '#endif'].join('\n');
    expect(preprocessGlsl(expression, new Map([['dN', '2']]))).toBe('e');
  });

  it('counts samplers and varyings', () => {
    const [, shaderCode] = GraphicsShaders[0];
    const values = {
      dRenderVariant: 'colorBlended',
      dGeometryType: 'mesh',
      dColorType: 'texture',
      dColorMarker: true,
      dLightCount: 1,
    };
    const defines = getProgramDefines(shaderCode, values, WebGL2PrefixDefines);
    defines.set('__VERSION__', '300');
    const vert = preprocessProgram(getGlslProgram(shaderCode.name).vert, defines, values);
    expect(getActiveSamplers(vert, getSchema(GraphicsShaders[0][2], false))).toEqual(
      expect.arrayContaining(['tColor', 'tMarker']),
    );
    expect(getVaryingTypes(vert)).toEqual(expect.arrayContaining(['vec4', 'vec3']));
    expect(getPackedVaryingRows(['mat4', 'vec4', 'vec3', 'vec3', 'float', 'float', 'float', 'vec2'])).toBe(8);
  });
  for (const [name, shaderCode, renderableSchema, geometryType] of GraphicsShaders) {
    const schema = getSchema(renderableSchema, true);
    for (const variant of RenderVariants) {
      const values: { [k: string]: string | boolean | number } = {
        dRenderVariant: variant,
        dGeometryType: geometryType,
      };
      for (const k of Object.keys(FeatureDefines)) {
        if (schema[k]?.type === 'define') values[k] = FeatureDefines[k];
      }
      const defines = getProgramDefines(shaderCode, values, WebGL2PrefixDefines);
      defines.set('__VERSION__', '300');
      const vert = preprocessProgram(getGlslProgram(shaderCode.name).vert, defines, values);
      const frag = preprocessProgram(getGlslProgram(shaderCode.name).frag, defines, values);
      const id = `${name} ${variant}`;

      it(`declared varyings: ${id}`, () => {
        // WebGL1 without extensions takes different branches (e.g. picking without draw buffers)
        const webgl1 = getProgramDefines(shaderCode, values, []);
        webgl1.set('__VERSION__', '100');
        for (const source of [vert, frag]) {
          expect(getUndeclaredVaryings(source)).toEqual([]);
        }
        expect(getUndeclaredVaryings(preprocessProgram(getGlslProgram(shaderCode.name).vert, webgl1, values))).toEqual(
          [],
        );
        expect(getUndeclaredVaryings(preprocessProgram(getGlslProgram(shaderCode.name).frag, webgl1, values))).toEqual(
          [],
        );
      });

      it(`textures per stage: ${id}`, () => {
        const limit = Math.min(WebGL2Minimum.texturesPerStage, WebGPUDefault.sampledTexturesPerStage);
        expect(getActiveSamplers(vert, schema).length).toBeLessThanOrEqual(limit);
        expect(getActiveSamplers(frag, schema).length).toBeLessThanOrEqual(limit);
      });

      it(`varyings: ${id}`, () => {
        const types = getVaryingTypes(vert);
        expect(getPackedVaryingRows(types)).toBeLessThanOrEqual(WebGL2Minimum.varyingVectors);
        const exceeds = getInterStageLocations(types) > WebGPUDefault.interStageVariables;
        expect(exceeds).toBe(KnownInterStageExceedances.has(id));
      });

      it(`webgpu binding layout: ${id}`, () => {
        const lookup = getDefineLookup(defines.keys());
        const layout = getBindingLayout(schema, lookup);
        const count = layout.global.textures.length + layout.material.textures.length + layout.object.textures.length;
        expect(count > WebGPUDefault.sampledTexturesPerStage).toBe(KnownBindingLayoutExceedances.has(id));
        // textures excluded by `usedIf` must not be used
        for (const k of [...getActiveSamplers(vert, schema), ...getActiveSamplers(frag, schema)]) {
          const spec = schema[k];
          if (spec?.type === 'texture') expect([k, isUsed(spec, lookup)]).toEqual([k, true]);
        }
      });
    }
  }
});

describe('minimum limits: compute shaders', () => {
  const textureLimit = Math.min(WebGL2Minimum.texturesPerStage, WebGPUDefault.sampledTexturesPerStage);

  for (const [name, shaderCode, schema] of ComputeShaders) {
    it(`attributes: ${name}`, () => {
      const layout = getAttributeLayout(schema);
      expect(layout.size).toBeLessThanOrEqual(Math.min(WebGL2Minimum.vertexAttribs, WebGPUDefault.vertexAttributes));
      expect(layout.entries.length).toBeLessThanOrEqual(WebGPUDefault.vertexBuffers);
    });

    it(`block sizes: ${name}`, () => {
      for (const block of UniformBlocks) {
        const size = getUniformBlockLayout(block, schema)?.size ?? 0;
        expect(size).toBeLessThanOrEqual(
          Math.min(WebGL2Minimum.uniformBlockSize, WebGPUDefault.uniformBufferBindingSize),
        );
      }
    });

    it(`per program: ${name}`, () => {
      const problems = new Set<string>();
      for (const values of getComputeDefineVariants(schema)) {
        for (const isWebGL2 of [true, false]) {
          const defines = getProgramDefines(shaderCode, values, getComputePrefixDefines(shaderCode, isWebGL2));
          defines.set('__VERSION__', isWebGL2 ? '300' : '100');
          const vert = preprocessProgram(getGlslProgram(shaderCode.name).vert, defines, values);
          const frag = preprocessProgram(getGlslProgram(shaderCode.name).frag, defines, values);
          const version = `webgl${isWebGL2 ? 2 : 1}`;
          if (!/void\s+main\s*\(/.test(vert) || !/void\s+main\s*\(/.test(frag)) {
            problems.add(`${version}: no main after preprocessing`);
          }

          for (const v of [...getUndeclaredVaryings(vert), ...getUndeclaredVaryings(frag)]) {
            problems.add(`${version}: undeclared varying '${v}'`);
          }
          const vertSamplers = getActiveSamplers(vert, schema);
          const fragSamplers = getActiveSamplers(frag, schema);
          if (Math.max(vertSamplers.length, fragSamplers.length) > textureLimit) {
            problems.add(`${version}: too many textures per stage`);
          }
          const types = getVaryingTypes(vert);
          if (getPackedVaryingRows(types) > WebGL2Minimum.varyingVectors) {
            problems.add(`${version}: too many varyings`);
          }
          if (getInterStageLocations(types) > WebGPUDefault.interStageVariables) {
            problems.add(`${version}: too many inter-stage variables`);
          }

          const lookup = getDefineLookup(getProgramDefines(shaderCode, values, []).keys());
          const layout = getBindingLayout(schema, lookup);
          const count = layout.global.textures.length + layout.material.textures.length + layout.object.textures.length;
          if (count > WebGPUDefault.sampledTexturesPerStage) problems.add('too many bound textures');
          for (const k of [...vertSamplers, ...fragSamplers]) {
            const spec = schema[k];
            if (spec?.type !== 'texture') problems.add(`${version}: sampler '${k}' not in schema`);
            else if (!isUsed(spec, lookup)) problems.add(`'${k}' used but excluded by usedIf`);
          }
        }
      }
      expect([...problems]).toEqual([]);
    });
  }
});
