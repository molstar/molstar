/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { BaseSchema, GlobalDefineSchema, GlobalTextureSchema, GlobalUniformSchema } from '../renderable/schema.js';
import { getGlslProgram } from '../webgl/programs.js';
import { getWgslProgram } from '../webgpu/shader/programs.js';
import { GraphicsShaders, resolveGeometryType } from './graphics-shaders.js';

// The GLSL and WGSL sources of each graphics program kept in step: the schema names (defines, uniforms, textures)
// each backend's source reads, per geometry, the same but for the differences listed here. Not a translator: a
// feature added to one backend's shaders and not the other's shows as a new name on one side. The WGSL source is
// taken with every feature of its variants on (its defines' lookup all true), the superset of what it can read.

/**
 * The names one backend reads and the other doesn't, by program, as intended:
 * - WebGL mechanics: `uDrawId` (multi-draw), `uPickType`, `uClipObjectCount` and `dLightCount` (loop bounds as uniforms
 *   or defines, WGSL loops over arrays), the texture sizes WebGL samples data textures by (`uGeoTexDim`,
 *   `uBrickTexDim`, `uGridTexDim`), `uInvView` and `uInvariantBoundingSphere` (WGSL derives them), `dClipping` of
 *   direct volumes (`#undef`, no clipping layers there on either backend);
 * - grids (`t*Grid`, `u*GridTexDim`): only meshes have them (color smoothing); the shared GLSL color chunk reads them
 *   for every geometry, WebGPU reports them unsupported for the others (`getUnsupported`);
 * - WGSL: direct volumes' `uInvView` (GLSL derives it), the sphere raster's uniforms (`uRaster*`, `uViewport`).
 */
const Grids = [
  'tColorGrid',
  'tEmissiveGrid',
  'tOverpaintGrid',
  'tSubstanceGrid',
  'tTransparencyGrid',
  'uColorGridTexDim',
  'uEmissiveGridTexDim',
  'uOverpaintGridTexDim',
  'uSubstanceGridTexDim',
  'uTransparencyGridTexDim',
];
const Expected: { [program: string]: { glsl: string[]; wgsl: string[] } } = {
  mesh: { glsl: ['dLightCount', 'uClipObjectCount', 'uDrawId', 'uPickType'], wgsl: [] },
  'texture-mesh': {
    glsl: ['dLightCount', 'uClipObjectCount', 'uDrawId', 'uGeoTexDim', 'uPickType'],
    wgsl: [],
  },
  spheres: {
    glsl: ['dLightCount', ...Grids, 'uClipObjectCount', 'uDrawId', 'uPickType'],
    wgsl: ['uRasterChunkBase', 'uRasterChunkCount', 'uRasterLevels', 'uRasterMaxRadius', 'uRasterWidth', 'uViewport'],
  },
  cylinders: {
    glsl: ['dLightCount', ...Grids, 'uClipObjectCount', 'uDrawId', 'uPickType'],
    wgsl: [],
  },
  points: { glsl: [...Grids, 'uClipObjectCount', 'uDrawId', 'uInvView', 'uPickType'], wgsl: [] },
  lines: { glsl: [...Grids, 'uClipObjectCount', 'uDrawId', 'uInvView', 'uPickType'], wgsl: [] },
  text: { glsl: [...Grids, 'uClipObjectCount', 'uDrawId', 'uInvView', 'uPickType'], wgsl: [] },
  image: {
    glsl: ['uClipObjectCount', 'uDrawId', 'uInvView', 'uInvariantBoundingSphere', 'uPickType'],
    wgsl: [],
  },
  'direct-volume': {
    glsl: ['dClipping', 'dLightCount', 'uBrickTexDim', 'uClipObjectCount', 'uGridTexDim', 'uPickType'],
    wgsl: ['uInvView'],
  },
};

/** The schema names `source` reads: whole words, a define's value suffix (`dX_value`) as the define; comments not */
function readNames(source: string, names: ReadonlySet<string>) {
  const code = source.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, '');
  const out = new Set<string>();
  for (const m of code.matchAll(/\b([dut][A-Z]\w*)\b/g)) {
    const name = m[1];
    if (names.has(name)) out.add(name);
    else {
      const define = name.replace(/_\w+$/, '');
      if (names.has(define)) out.add(define);
    }
  }
  return out;
}

describe('GLSL and WGSL in step', () => {
  for (const [name, shaderCode, schema, geometryType] of GraphicsShaders) {
    it(`${name}: read the same schema names, but for the listed differences`, () => {
      const names = new Set([
        ...Object.keys(schema),
        ...Object.keys(BaseSchema),
        ...Object.keys(GlobalUniformSchema),
        ...Object.keys(GlobalTextureSchema),
        ...Object.keys(GlobalDefineSchema),
      ]);
      const glsl = getGlslProgram(shaderCode.name);
      const wgsl = getWgslProgram(shaderCode.name);
      expect(wgsl, `no WGSL program '${shaderCode.name}'`).toBeDefined();
      const wgslSource = typeof wgsl!.source === 'string' ? wgsl!.source : wgsl!.source(() => true);
      const fromGlsl = readNames(
        resolveGeometryType(glsl.vert, geometryType) + resolveGeometryType(glsl.frag, geometryType),
        names,
      );
      const fromWgsl = readNames(resolveGeometryType(wgslSource, geometryType), names);
      const onlyGlsl = [...fromGlsl].filter((n) => !fromWgsl.has(n)).sort();
      const onlyWgsl = [...fromWgsl].filter((n) => !fromGlsl.has(n)).sort();
      // a new difference fails, and so does a listed one gone (the list kept to what is so)
      expect({ glsl: onlyGlsl, wgsl: onlyWgsl }).toEqual({
        glsl: [...Expected[name].glsl].sort(),
        wgsl: [...Expected[name].wgsl].sort(),
      });
    });
  }
});
