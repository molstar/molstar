/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// What an extension uses of the WebGL backend: its programs' GLSL registered (beside their WGSL, `gl/extend/webgpu`).
// The backend folder (`gl/webgl/`) is internal to the package (its exports source-only): its GLSL sources, chunks and
// assembly and its WebGL resources stay internal, an extension brings its own GLSL

export { type FragOutTypes, type GlslProgram, hasGlslProgram, registerGlslProgram } from '../webgl/program-registry.js';
