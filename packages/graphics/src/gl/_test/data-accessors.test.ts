/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { getGlslProgram } from '../webgl/programs.js';
import { getWgslProgram, resolveWgslSource } from '../webgpu/shader/programs.js';
import { GraphicsShaders } from './graphics-shaders.js';

/** `source` without comments: they may name the arrays */
function stripComments(source: string) {
  return source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, '');
}

/** The names of `names` that `source` references other than through their accessors (`read_<name>`) */
function getDirectReads(source: string, names: readonly string[]) {
  const code = stripComments(source);
  return names.filter((name) => new RegExp(`(?<!read_)\\b${name}\\b`).test(code));
}

// The data arrays (`DataSpec`) have no layout of their own: the backends lay them out (WebGL data textures, WebGPU
// storage pages), so shaders read them only through the generated accessors
describe('data arrays', () => {
  for (const [name, shaderCode, schema] of GraphicsShaders) {
    const names = Object.keys(schema).filter((k) => schema[k].type === 'data');
    it(`are read through their accessors: ${name}`, () => {
      expect(names.length).toBeGreaterThan(0);
      const glsl = getGlslProgram(shaderCode.name);
      expect(getDirectReads(glsl.vert, names), 'GLSL vertex').toEqual([]);
      expect(getDirectReads(glsl.frag, names), 'GLSL fragment').toEqual([]);
      const wgsl = getWgslProgram(shaderCode.name);
      if (wgsl) expect(getDirectReads(resolveWgslSource(wgsl.source, {}), names), 'WGSL').toEqual([]);
    });
  }
});
