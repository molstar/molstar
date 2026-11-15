/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

// The WGSL compute kernels by name, without the package's own (`kernels.ts` registers those): what an extension
// registers into (`gl/extend/webgpu`) imports none of the built-in sources

const Kernels = new Map<string, string>();

/**
 * Registers the WGSL source of the kernel `name`, its bindings declared from its schema (`getWgslKernelDeclarations`);
 * throws where another source has the name
 */
export function registerWgslKernel(name: string, source: string) {
  if (Kernels.has(name) && Kernels.get(name) !== source) throw new Error(`kernel '${name}' is already registered`);
  Kernels.set(name, source);
}

export function getWgslKernel(name: string): string | undefined {
  return Kernels.get(name);
}

export function getWgslKernelNames(): string[] {
  return [...Kernels.keys()];
}
