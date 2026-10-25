/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import {
  GlobalDefineSchema,
  GlobalTextureSchema,
  GlobalUniformSchema,
  getDefineLookup,
  InternalSchema,
  type RenderableSchema,
} from '../renderable/schema.js';
import type { ShaderCode } from '../shader-code.js';
import { preprocessShader, ShaderPreprocessorError, unrollLoops } from '../shader-preprocessor.js';
import { getShaderDeclarations } from '../webgl/shader-declarations.js';
import { ComputeShaders, getComputeDefineVariants, getComputePrefixDefines } from './compute-shaders.js';
import { GraphicsShaders, getProgramDefines, preprocessGlsl } from './graphics-shaders.js';
import { getGlslProgram } from '../webgl/programs.js';

const lines = (...l: string[]) => l.join('\n');
const run = (source: string, defines: [string, string][] = [], expandMacros = false) =>
  preprocessShader(source, new Map(defines), { expandMacros });

describe('shader preprocessor', () => {
  it('conditionals', () => {
    const source = lines(
      '#if defined(dA) && !defined(dB)',
      'a',
      '#elif dN > 1',
      'b',
      '#else',
      'c',
      '#endif',
      '#ifdef dB',
      '  #ifndef dA',
      '  d',
      '  #endif',
      '#endif',
    );
    expect(run(source, [['dA', '']])).toBe('a');
    expect(run(source, [['dN', '2']])).toBe('b');
    expect(run(source, [['dN', '1']])).toBe('c');
    expect(
      run(source, [
        ['dN', '0'],
        ['dB', ''],
      ]),
    ).toBe(lines('c', '  d'));
    expect(
      run(lines('#if defined dA || (dN == 2 && dM != 3)', 'x', '#endif'), [
        ['dN', '2'],
        ['dM', '4'],
      ]),
    ).toBe('x');
    expect(run(lines('#if dN >= 2 && dN <= 2 && !(dN < 2)', 'x', '#endif'), [['dN', '2']])).toBe('x');
  });

  it('skips inactive branches without evaluating them', () => {
    // undefined names in skipped branches are fine, like in the C preprocessor
    expect(run(lines('#if 0', '#if dUndefined > 1', 'x', '#endif', '#elif 1', 'y', '#endif'))).toBe('y');
    expect(run(lines('#if 1', 'a', '#elif dUndefined', 'b', '#endif'))).toBe('a');
  });

  it('keeps or expands macros', () => {
    const source = lines(
      '#define saturate(a) clamp(a, 0.0, 1.0)',
      '#define SCALE 2.0',
      'float x = saturate(v * SCALE) + dN;',
      '#undef SCALE',
      'float y = SCALE;',
    );
    expect(run(source, [['dN', '3']])).toBe(source);
    expect(run(source, [['dN', '3']], true)).toBe(lines('float x = clamp(v * 2.0, 0.0, 1.0) + 3;', 'float y = SCALE;'));
    // nested invocations, parentheses in arguments, no recursion, numbers are not names
    const nested = lines(
      '#define F(a, b) (a + b)',
      '#define G(a) F(a, a)',
      '#define R R + 1',
      'x = G(f(1, 2)) * R + 1e5 + F;',
    );
    expect(run(nested, [], true)).toBe('x = (f(1, 2) + f(1, 2)) * R + 1 + 1e5 + F;');
    // inactive defines don't apply
    expect(run(lines('#ifdef dA', '#define X 1', '#endif', 'X'), [], true)).toBe('X');
  });

  it('unrolls loops', () => {
    const source = lines(
      '#pragma unroll_loop_start',
      'for (int i = 0; i < 2; ++i) {',
      '  a[i] += UNROLLED_LOOP_INDEX;',
      '}',
      '#pragma unroll_loop_end',
    );
    const code = (s: string) =>
      s
        .split('\n')
        .map((l) => l.trim())
        .filter((l) => l);
    expect(code(unrollLoops(source))).toEqual(['a[0] += 0;', 'a[1] += 1;']);
    expect(code(run(source))).toEqual(['a[0] += 0;', 'a[1] += 1;']);
  });

  it('is strict', () => {
    const fails = (source: string, message: RegExp, defines: [string, string][] = []) => {
      expect(() => run(source, defines)).toThrow(ShaderPreprocessorError);
      expect(() => run(source, defines)).toThrow(message);
    };
    fails(lines('#if dA', '#endif'), /undefined name 'dA'/);
    fails(lines('#if dA', '#endif'), /'dA' is not an integer/, [['dA', '1.5']]);
    fails(lines('#if dN + 1 > 2', '#endif'), /unsupported token/, [['dN', '1']]);
    fails(lines('#if 1 ?', '#endif'), /unsupported token/);
    fails(lines('#if (1', '#endif'), /unexpected end/);
    fails(lines('#if 1 1', '#endif'), /unexpected '1'/);
    fails(lines('#if', '#endif'), /empty expression/);
    fails(lines('#ifdef', '#endif'), /expects a name/);
    fails('#if 1', /missing #endif/);
    fails('#endif', /#endif without #if/);
    fails(lines('#if 1', '#else', '#else', '#endif'), /#else after #else/);
    fails(lines('#if 1', '#else', '#elif 1', '#endif'), /#elif after #else/);
    fails('#version 300 es', /unsupported directive '#version'/);
    fails('#extension GL_EXT_frag_depth : enable', /unsupported directive/);
    fails('#pragma optimize(on)', /unsupported pragma/);
    fails('#define STR(x) #x', /'#' and '##' are not supported/);
    expect(() => run(lines('#define F(a, b) a', 'F(1)'), [], true)).toThrow(/expects 2 arguments, got 1/);
    expect(() => run(lines('#define F(a) a', 'F(1'), [], true)).toThrow(/unterminated invocation/);
  });
});

//

const RenderVariants = ['colorBlended', 'colorWboit', 'colorDpoit', 'pick', 'depth', 'marking', 'emissive', 'tracing'];
/** all define variants for these, the base variant for the others */
const FullRenderVariants = new Set(['colorBlended', 'pick']);

type DefineValues = { [k: string]: string | boolean | number };

/** the count replacement of `shader-code.ts` (`replaceCounts`), applied before preprocessing */
function replaceCounts(source: string, values: DefineValues) {
  for (const name of ['dLightCount']) {
    if (values[name] !== undefined) source = source.replace(new RegExp(name, 'g'), `${values[name]}`);
  }
  return source;
}

function withoutDefines(source: string) {
  return source
    .split('\n')
    .filter((l) => !/^\s*#\s*(define|undef)\b/.test(l))
    .join('\n');
}

/** checks a program variant against the lenient reference preprocessor of the shader tests */
function checkVariant(shaderCode: ShaderCode, schema: RenderableSchema, values: DefineValues, isWebGL2: boolean) {
  const defines = getProgramDefines(shaderCode, values, getComputePrefixDefines(shaderCode, isWebGL2));
  defines.set('__VERSION__', isWebGL2 ? '300' : '100');
  const declarations = getShaderDeclarations(schema, isWebGL2, getDefineLookup(defines.keys()));
  const sources = [
    declarations.vert,
    getGlslProgram(shaderCode.name).vert,
    declarations.frag,
    getGlslProgram(shaderCode.name).frag,
  ].map((s) => unrollLoops(replaceCounts(s, values)));
  for (const source of sources) {
    const kept = preprocessShader(source, defines);
    expect(withoutDefines(kept)).toBe(withoutDefines(preprocessGlsl(source, defines)));
    const expanded = preprocessShader(source, defines, { expandMacros: true });
    expect(expanded).not.toMatch(/^\s*#/m);
    // number defines are substituted
    for (const [k, v] of defines) {
      if (/^\d+$/.test(v)) expect([k, new RegExp(`\\b${k}\\b`).test(expanded)]).toEqual([k, false]);
    }
  }
}

// every shader in every variant: the mesh is the quick lane's representative, the others are medium (AGENTS.md, Tests)
describe('shader preprocessor: inventory', () => {
  for (const [name, shaderCode, renderableSchema, geometryType] of GraphicsShaders) {
    const schema: RenderableSchema = {
      ...GlobalUniformSchema,
      ...GlobalTextureSchema,
      ...GlobalDefineSchema,
      ...InternalSchema,
      ...renderableSchema,
    };
    const defineVariants = getComputeDefineVariants(schema);
    for (const renderVariant of RenderVariants) {
      it(`${name} ${renderVariant}`, name === 'mesh' ? {} : { tags: ['medium'] }, () => {
        const variants = FullRenderVariants.has(renderVariant) ? defineVariants : defineVariants.slice(0, 1);
        for (const variant of variants) {
          const values = { ...variant, dRenderVariant: renderVariant, dGeometryType: geometryType };
          for (const isWebGL2 of [true, false]) checkVariant(shaderCode, schema, values, isWebGL2);
        }
      });
    }
  }

  for (const [name, shaderCode, schema] of ComputeShaders) {
    it(`compute ${name}`, () => {
      for (const values of getComputeDefineVariants(schema)) {
        for (const isWebGL2 of [true, false]) checkVariant(shaderCode, schema, values, isWebGL2);
      }
    });
  }
});
