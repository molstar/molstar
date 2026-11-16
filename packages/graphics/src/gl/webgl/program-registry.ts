/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The GLSL programs by name, without the package's own (`programs.ts` registers those): what an extension registers
// into (`gl/extend/webgl`) imports none of the built-in sources

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
