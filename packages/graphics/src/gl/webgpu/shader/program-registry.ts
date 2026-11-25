/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The WGSL programs by name, without the package's own (`programs.ts` registers those): what an extension registers
// into (`gl/extend/webgpu`) imports none of the built-in sources

import type { ProgramVariant } from '../../render-item.js';
import type { DefineLookup } from '../../renderable/schema.js';
import type { DefineValues } from '../../shader-code.js';

/**
 * A WGSL program, by the name of the `ShaderCode` it replaces. `source` may depend on the variant's defines (a lookup
 * of the flattened names, `dX` and `dX_value`): what a variant turns off (a feature's helpers, the reads of its
 * textures) is left out of the source, so that it is not bound; defines inside the source are for small branches
 */
export interface WgslProgram {
  readonly source: string | ((defines: DefineLookup) => string);
  /** render variants the source implements */
  readonly variants: readonly ProgramVariant[];
  /** the first feature selected by `defines` that the source doesn't implement */
  readonly getUnsupported?: (defines: DefineValues) => string | undefined;
}

const Programs = new Map<string, WgslProgram>();
/** the package's own programs, which `unregisterWgslProgram` keeps */
const BuiltIn = new Set<string>();

/** The package's own programs (`programs.ts`), registered once */
export function registerBuiltInWgslPrograms(programs: { readonly [name: string]: WgslProgram }) {
  for (const name of Object.keys(programs)) {
    Programs.set(name, programs[name]);
    BuiltIn.add(name);
  }
}

/**
 * Registers the WGSL program of the `ShaderCode` named `name` (an extension's program, beside its GLSL one,
 * `registerGlslProgram`): its draws on WebGPU then use it; throws where another program has the name
 */
export function registerWgslProgram(name: string, program: WgslProgram) {
  const registered = Programs.get(name);
  if (registered && registered !== program) throw new Error(`WGSL program '${name}' is already registered`);
  Programs.set(name, program);
}

/** Removes a program `registerWgslProgram` added (a test restoring the registry); the package's own stay */
export function unregisterWgslProgram(name: string) {
  if (!BuiltIn.has(name)) Programs.delete(name);
}

export function getWgslProgram(name: string): WgslProgram | undefined {
  return Programs.get(name);
}

export function hasWgslProgram(name: string) {
  return Programs.has(name);
}

/** The names of the registered programs (e.g. for a check that every one compiles) */
export function getWgslProgramNames(): string[] {
  return [...Programs.keys()];
}
