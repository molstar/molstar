/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { afterEach, describe, expect, it } from 'vitest';
import type { ShaderCode } from '../shader-code.js';
import { hasWgslProgram, registerWgslProgram, type WgslProgram } from '../extend/webgpu.js';
import { createWebGPUContext } from '../webgpu/context.js';
import { getWgslProgram, unregisterWgslProgram } from '../webgpu/shader/programs.js';
import { createStubDevice } from './stub-gpu.js';

const Name = 'test-extension-program';
const program: WgslProgram = { source: '// an extension program', variants: ['color'] };

describe('registering a WGSL program', () => {
  // the registry is module state other files share
  afterEach(() => unregisterWgslProgram(Name));

  it('makes a program of the name drawn on WebGPU', () => {
    const ctx = createWebGPUContext(createStubDevice(), { width: 16, height: 16 });
    const shaderCode = { name: Name } as ShaderCode;
    expect(ctx.isProgramSupported(shaderCode)).toBe(false);
    registerWgslProgram(Name, program);
    expect(hasWgslProgram(Name)).toBe(true);
    expect(getWgslProgram(Name)).toBe(program);
    expect(ctx.isProgramSupported(shaderCode)).toBe(true);
    ctx.destroy();
  });

  it('keeps a name to its first program', () => {
    registerWgslProgram(Name, program);
    // the same program again is fine (a module evaluated twice), another one is not
    expect(() => registerWgslProgram(Name, program)).not.toThrow();
    expect(() => registerWgslProgram(Name, { ...program })).toThrow(/already registered/);
    expect(() => registerWgslProgram('mesh', program)).toThrow(/already registered/);
  });

  it('removes only what was registered, not the package programs', () => {
    unregisterWgslProgram('mesh');
    expect(hasWgslProgram('mesh')).toBe(true);
  });
});
