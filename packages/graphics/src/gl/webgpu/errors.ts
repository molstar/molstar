/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/** A feature the WebGPU backend doesn't implement yet */
export class WebGPUUnsupportedError extends Error {
  constructor(feature: string) {
    super(`${feature} is not supported by the WebGPU backend yet`);
    this.name = 'WebGPUUnsupportedError';
  }
}

export function notSupported(feature: string): never {
  throw new WebGPUUnsupportedError(feature);
}
