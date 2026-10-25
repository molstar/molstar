/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { getGrid3dComputeFrag, getGrid3dComputeSchema } from '../webgl/compute/grid3d.js';
import {
  AttributeSpec,
  BaseSchema,
  getDefineLookup,
  GlobalDefineSchema,
  GlobalTextureSchema,
  GlobalUniformSchema,
  InternalSchema,
  type RenderableSchema,
  SegmentSchema,
  TextureSpec,
  UniformSpec,
} from '../renderable/schema.js';
import { ShaderCode } from '../shader-code.js';
import {
  stripAttributeDeclarations,
  stripSamplerDeclarations,
  stripUniformDeclarations,
} from '../webgl/shader-code.js';
import { quad_vert } from '../webgl/shader/quad.vert.js';
import { getUniformBlockLayout, UniformBlockInfo, UniformBlocks } from '../uniform-layout.js';
import { getShaderDeclarations } from '../webgl/shader-declarations.js';
import {
  getWgslDrawBinding,
  getWgslDrawLayout,
  getWgslUniformDeclarations,
  WgslUniformBlockInstance,
} from '../webgpu/uniform-declarations.js';
import { ComputeShaders, getComputeDefineVariants } from './compute-shaders.js';
import { GraphicsShaders, getProgramDefines } from './graphics-shaders.js';
import { getGlslProgram, registerGlslProgram } from '../webgl/programs.js';
import { getInstanceDataSchema } from '../webgl/render-item.js';

/**
 * Removes inactive branches of conditionals that only test `dGeometryType_*`,
 * all branches of other conditionals are kept.
 */
function resolveGeometryType(source: string, geometryType: string) {
  const evaluate = (condition: string): boolean | undefined => {
    const m = condition
      .trim()
      .match(/^(!?)\s*(?:defined\s*\(\s*(dGeometryType_\w+)\s*\)|defined\s+(dGeometryType_\w+))$/);
    if (!m) return undefined;
    return ((m[2] ?? m[3]) === `dGeometryType_${geometryType}`) !== (m[1] === '!');
  };
  // `active` of the current branch, `taken` if any branch of a resolved conditional was active
  const stack: { resolved: boolean; active: boolean; taken: boolean }[] = [];
  const lines: string[] = [];
  for (const line of source.split('\n')) {
    const d = line.trim().match(/^#\s*(ifdef|ifndef|if|elif|else|endif)\b(.*)$/);
    if (!d) {
      if (stack.every((f) => f.active)) lines.push(line);
      continue;
    }
    const [, directive, rest] = d;
    const top = stack[stack.length - 1];
    if (directive === 'ifdef' || directive === 'ifndef' || directive === 'if') {
      const condition = directive === 'if' ? rest : `${directive === 'ifndef' ? '!' : ''}defined(${rest.trim()})`;
      const value = evaluate(condition);
      stack.push(
        value === undefined
          ? { resolved: false, active: true, taken: false }
          : { resolved: true, active: value, taken: value },
      );
    } else if (directive === 'elif') {
      if (top.resolved) {
        const value = evaluate(rest) ?? true;
        top.active = !top.taken && value;
        top.taken ||= top.active;
      }
    } else if (directive === 'else') {
      if (top.resolved) top.active = !top.taken;
    } else {
      stack.pop();
    }
  }
  return lines.join('\n');
}

function getIdentifiers(source: string, prefix: 'u' | 't' | 'a') {
  const code = source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, '');
  return new Set(code.match(new RegExp(`\\b${prefix}[A-Z]\\w*`, 'g')) ?? []);
}

/** referenced in a shared chunk only in `bumpEnabled` branches, defined by shaders of geometries with these uniforms */
const ConditionalUniforms = new Set(['uBumpFrequency', 'uBumpAmplitude']);

const Guards: [string, 'u' | 't' | 'a', RenderableSchema['k']['type'], Set<string>][] = [
  ['uniforms', 'u', 'uniform', ConditionalUniforms],
  ['textures', 't', 'texture', new Set()],
  // `aColorMode` is cylinders-only (`dDualColor`)
  ['attributes', 'a', 'attribute', new Set(['aColorMode'])],
];

describe('shader identifiers', () => {
  // all identifiers referenced in any preprocessor branch of the shader code must be in the schema,
  // uniform declarations are generated from it (textures and attributes will be too)
  for (const [name, shaderCode, renderableSchema, geometryType] of GraphicsShaders) {
    // with what the WebGL backend adds for merged render items (`dInstanceData`)
    const schema = getInstanceDataSchema({
      ...GlobalUniformSchema,
      ...GlobalTextureSchema,
      ...GlobalDefineSchema,
      ...InternalSchema,
      ...renderableSchema,
      ...SegmentSchema,
    });
    const vert = resolveGeometryType(getGlslProgram(shaderCode.name).vert, geometryType);
    const frag = resolveGeometryType(getGlslProgram(shaderCode.name).frag, geometryType);
    for (const [kind, prefix, type, allowed] of Guards) {
      it(`${kind}: ${name}`, () => {
        const used = new Set([...getIdentifiers(vert, prefix), ...getIdentifiers(frag, prefix)]);
        const missing = [...used].filter((k) => schema[k]?.type !== type && !allowed.has(k));
        expect(missing).toEqual([]);
      });
    }
  }

  for (const [name, shaderCode, schema] of ComputeShaders) {
    for (const [kind, prefix, type] of Guards) {
      it(`${kind}: ${name}`, () => {
        const used = [
          ...getIdentifiers(getGlslProgram(shaderCode.name).vert, prefix),
          ...getIdentifiers(getGlslProgram(shaderCode.name).frag, prefix),
        ];
        const missing = [...new Set(used)].filter((k) => schema[k]?.type !== type);
        expect(missing).toEqual([]);
      });
    }
  }
});

describe('generated declarations', () => {
  // shaders must not hand-declare samplers, attributes and uniforms that are generated from the schema
  const shaders: [string, ShaderCode, RenderableSchema][] = [
    ...GraphicsShaders.map(
      ([name, shaderCode, schema]) =>
        [
          name,
          shaderCode,
          {
            ...GlobalUniformSchema,
            ...GlobalTextureSchema,
            ...GlobalDefineSchema,
            ...InternalSchema,
            ...schema,
            ...SegmentSchema,
          },
        ] as [string, ShaderCode, RenderableSchema],
    ),
    ...ComputeShaders,
  ];
  for (const [name, shaderCode, schema] of shaders) {
    it(`no hand-written declarations: ${name}`, () => {
      const source = `${getGlslProgram(shaderCode.name).vert}\n${getGlslProgram(shaderCode.name).frag}`;
      const declared = [
        ...source.matchAll(/uniform\s+(?:(?:lowp|mediump|highp)\s+)?\w*sampler\w+\s+(\w+)/g),
        ...source.matchAll(/^\s*attribute\s+(?:(?:lowp|mediump|highp)\s+)?\w+\s+(\w+)/gm),
        ...source.matchAll(/^\s*uniform\s+(?:(?:lowp|mediump|highp)\s+)?\w+\s+(\w+)\s*(?:\[[^\]]*\])?\s*;/gm),
      ].map((m) => m[1]);
      expect(declared.filter((k) => schema[k] !== undefined)).toEqual([]);
    });

    it(`array length defines: ${name}`, () => {
      const invalid = Object.entries(schema).filter(([, spec]) => {
        if (spec.type !== 'uniform' || spec.lengthDefine === undefined) return false;
        // a number define, or a loop-bound integer uniform (a define on WebGL1, the maximum `length` elsewhere)
        const define = schema[spec.lengthDefine];
        if (define?.type === 'uniform') return define.kind !== 'i' || !define.loopBound || spec.length === undefined;
        return define?.type !== 'define' || define.kind !== 'number';
      });
      expect(invalid.map(([k]) => k)).toEqual([]);
    });
  }

  it('sampler types', () => {
    const schema = GraphicsShaders.find(([name]) => name === 'direct-volume')![2];
    const lookup = getDefineLookup(['dGridTexType_3d']);
    expect(getShaderDeclarations(schema, true, lookup).frag).toContain('uniform highp sampler3D tGridTex;');
    expect(getShaderDeclarations(schema, false, getDefineLookup(['dGridTexType_2d'])).frag).toContain(
      'uniform highp sampler2D tGridTex;',
    );
    // no 3D textures on WebGL1, the hand-written declaration (in an inactive branch) is kept
    expect(getShaderDeclarations(schema, false, lookup).samplerTypes.has('tGridTex')).toBe(false);

    const integer = { tInt: TextureSpec('texture', 'rgba', 'float', 'nearest', undefined, { sampleType: 'sint' }) };
    expect(getShaderDeclarations(integer, true, lookup).frag).toBe('uniform highp isampler2D tInt;');
    expect(getShaderDeclarations(integer, false, lookup).frag).toBe('uniform highp sampler2D tInt;');
    expect(getShaderDeclarations(integer, true, lookup).samplerTypes.get('tInt')).toEqual(['isampler2D', 'sampler2D']);
  });

  it('attribute declarations', () => {
    const schema = {
      aPosition: AttributeSpec('float32', 3, 0),
      aTransform: AttributeSpec('float32', 16, 1),
      aGroup: AttributeSpec('float32', 1, 0),
    };
    const lookup = getDefineLookup([]);
    expect(getShaderDeclarations(schema, true, lookup).vert.split('\n')).toEqual([
      'layout(location = 0) in vec3 aPosition;',
      'layout(location = 1) in mat4 aTransform;',
      'layout(location = 5) in float aGroup;',
    ]);
    expect(getShaderDeclarations(schema, false, lookup).vert.split('\n')).toEqual([
      'attribute vec3 aPosition;',
      'attribute mat4 aTransform;',
      'attribute float aGroup;',
    ]);
    expect(getShaderDeclarations(schema, true, lookup).frag).toBe('');
  });

  it('strip', () => {
    const samplers = new Map([['tA', ['sampler2D', 'sampler3D']]]);
    const source = 'uniform highp sampler2D tA;\nuniform sampler3D tA;\nuniform sampler2D tB;';
    expect(stripSamplerDeclarations(source, samplers)).toBe('\n\nuniform sampler2D tB;');
    expect(() => stripSamplerDeclarations('uniform samplerCube tA;', samplers)).toThrow();

    const attributes = new Map([['aA', 'vec3']]);
    expect(stripAttributeDeclarations('attribute vec3 aA;\nattribute float aB;', attributes)).toBe(
      '\nattribute float aB;',
    );
    expect(() => stripAttributeDeclarations('attribute vec2 aA;', attributes)).toThrow();
  });

  // hand-written declarations of compute shaders are still stripped, their types must match the spec
  for (const [name, shaderCode, schema] of ComputeShaders) {
    it(`compute declarations match the schema: ${name}`, () => {
      const errors = new Set<string>();
      for (const values of getComputeDefineVariants(schema)) {
        const lookup = getDefineLookup(getProgramDefines(shaderCode, values, []).keys());
        for (const isWebGL2 of [true, false]) {
          const { uniformTypes, samplerTypes, attributeTypes } = getShaderDeclarations(schema, isWebGL2, lookup);
          for (const source of [getGlslProgram(shaderCode.name).vert, getGlslProgram(shaderCode.name).frag]) {
            try {
              stripUniformDeclarations(
                stripAttributeDeclarations(stripSamplerDeclarations(source, samplerTypes), attributeTypes),
                uniformTypes,
              );
            } catch (e) {
              errors.add(`webgl${isWebGL2 ? 2 : 1}: ${(e as Error).message}`);
            }
          }
        }
      }
      expect([...errors]).toEqual([]);
    });
  }
});

describe('grid3d compute', () => {
  const spec = {
    schema: {
      tData: TextureSpec('image-float32', 'rgba', 'float', 'nearest'),
      uN: UniformSpec('i'),
      uScale: UniformSpec('f'),
    },
    loopBounds: ['uN' as const],
    mainCode: 'float v = 0.0; for (int i = 0; i < uN; ++i) { v += texture2D(tData, vec2(float(i), 0.5)).r; }',
    returnCode: 'current + uScale * v',
    values: () => ({}) as any,
    cumulative: { states: () => [], update: () => {} },
  };
  const shaderCode = ShaderCode('grid3d-test');
  registerGlslProgram('grid3d-test', quad_vert, getGrid3dComputeFrag(spec));

  for (const isWebGL2 of [true, false]) {
    const schema = getGrid3dComputeSchema(spec, isWebGL2);
    const version = `webgl${isWebGL2 ? 2 : 1}`;

    it(`identifiers: ${version}`, () => {
      for (const [, prefix, type] of Guards) {
        const used = [
          ...getIdentifiers(getGlslProgram(shaderCode.name).vert, prefix),
          ...getIdentifiers(getGlslProgram(shaderCode.name).frag, prefix),
        ];
        // loop bounds are number defines on WebGL1
        const missing = [...new Set(used)].filter((k) => schema[k]?.type !== type && schema[k]?.type !== 'define');
        expect(missing).toEqual([]);
      }
    });

    it(`declarations: ${version}`, () => {
      const source = `${getGlslProgram(shaderCode.name).vert}\n${getGlslProgram(shaderCode.name).frag}`;
      expect(source).not.toMatch(/^\s*(uniform|attribute)\b/m);
      const lookup = getDefineLookup(isWebGL2 ? ['CUMULATIVE'] : ['CUMULATIVE', 'WEBGL1']);
      const { frag, samplerTypes } = getShaderDeclarations(schema, isWebGL2, lookup);
      expect([...samplerTypes.keys()].sort()).toEqual(['tCumulativeSum', 'tData']);
      expect(/\buN\b/.test(frag)).toBe(isWebGL2);
      expect(frag).toContain('uScale');
    });
  }
});

/** [align, size] of a WGSL type in the uniform address space */
function getWgslAlignSize(type: string): [number, number] {
  const array = type.match(/^array<(.+), (\d+)>$/);
  if (array) {
    const [align, size] = getWgslAlignSize(array[1]);
    const stride = Math.ceil(size / align) * align;
    if (stride % 16 !== 0) throw new Error(`uniform array stride must be a multiple of 16: ${type}`);
    return [Math.max(align, 16), stride * Number(array[2])];
  }
  if (/^(f32|i32|u32)$/.test(type)) return [4, 4];
  if (type.startsWith('vec2<')) return [8, 8];
  if (type.startsWith('vec3<')) return [16, 12];
  if (type.startsWith('vec4<')) return [16, 16];
  if (type === 'mat3x3<f32>') return [16, 48];
  if (type === 'mat4x4<f32>') return [16, 64];
  throw new Error(`unknown WGSL type '${type}'`);
}

/** member offsets according to the WGSL memory layout rules */
function getWgslOffsets(code: string, struct: string) {
  const body = code.match(new RegExp(`struct ${struct} \\{\\n([\\s\\S]*?)\\n\\};`))![1];
  let offset = 0;
  return body.split('\n').map((line) => {
    const [, name, type] = line.trim().match(/^(\w+): (.+),$/)!;
    const [align, size] = getWgslAlignSize(type);
    offset = Math.ceil(offset / align) * align;
    const member: [string, number] = [name, offset];
    offset += size;
    return member;
  });
}

describe('wgsl uniforms', () => {
  const ArraySchema: RenderableSchema = {
    s: UniformSpec('f'),
    arrF: UniformSpec('f[]', undefined, undefined, 3),
    arrV2: UniformSpec('v2[]', undefined, undefined, 2),
    arrV3: UniformSpec('v3[]', undefined, undefined, 2),
    arrM3: UniformSpec('m3[]', undefined, undefined, 2),
    after: UniformSpec('f'),
    arrB: UniformSpec('b[]', undefined, undefined, 2),
    uD: UniformSpec('i', undefined, 'draw'),
    uC: UniformSpec('f[]'),
  };
  const schemas: [string, RenderableSchema][] = [
    ['arrays', ArraySchema],
    ...GraphicsShaders.map(
      ([name, , schema]) =>
        [name, { ...GlobalUniformSchema, ...InternalSchema, ...schema }] as [string, RenderableSchema],
    ),
  ];

  for (const [name, schema] of schemas) {
    it(`layout matches std140: ${name}`, () => {
      const { code } = getWgslUniformDeclarations(schema);
      for (const block of UniformBlocks) {
        const layout = getUniformBlockLayout(block, schema);
        if (!layout) continue;
        const { name: struct, bindingIndex } = UniformBlockInfo[block];
        expect(getWgslOffsets(code, struct)).toEqual(layout.members.map((m) => [m.name, m.offset]));
        expect(code).toContain(
          `@group(${bindingIndex}) @binding(0) var<uniform> ${WgslUniformBlockInstance[block]}: ${struct};`,
        );
      }
    });
  }

  it('unsupported', () => {
    expect(getWgslUniformDeclarations(ArraySchema).unsupported).toEqual(['uC']);
    expect(getWgslUniformDeclarations({ ...GlobalUniformSchema, ...BaseSchema, ...SegmentSchema }).unsupported).toEqual(
      [],
    );
  });

  it('draw uniforms', () => {
    const schema = { ...GlobalUniformSchema, ...BaseSchema, ...SegmentSchema };
    const { code } = getWgslUniformDeclarations(schema);
    const layout = getWgslDrawLayout(schema)!;
    expect(layout.members.map((m) => m.name)).toEqual([
      'uDrawId',
      'uSolidInteriorPass',
      'uSolidInteriorClip',
      'uLod',
      'uLodLevel',
    ]);
    expect(getWgslOffsets(code, 'Draw')).toEqual(layout.members.map((m) => [m.name, m.offset]));
    expect(code).toContain(`@group(2) @binding(${getWgslDrawBinding(schema)}) var<uniform> uDraw: Draw;`);
  });
});
