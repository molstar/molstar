/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// What an extension uses of the WebGPU backend: its programs' WGSL (beside their GLSL, `gl/extend/webgl`) and its
// compute kernels registered. The backend folder (`gl/webgpu/`) is internal to the package (its exports source-only)

export { hasWgslProgram, registerWgslProgram, type WgslProgram } from '../webgpu/shader/program-registry.js';
export { registerWgslKernel } from '../webgpu/shader/kernel-registry.js';
