/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it } from 'vitest';
import { getAdapterMaxLimits } from '../webgpu/context.js';

/** Supported limits as browsers expose them: getters of the prototype */
class Limits {
  get maxTextureDimension2D() {
    return 16384;
  }
  get maxStorageBufferBindingSize() {
    return 2147483644;
  }
  get maxComputeWorkgroupsPerDimension() {
    return 65535;
  }
  get minUniformBufferOffsetAlignment() {
    return 64;
  }
  get minStorageBufferOffsetAlignment() {
    return 32;
  }
}
// enumerable as WebIDL attributes are
for (const key of Object.getOwnPropertyNames(Limits.prototype)) {
  if (key === 'constructor') continue;
  Object.defineProperty(Limits.prototype, key, {
    ...Object.getOwnPropertyDescriptor(Limits.prototype, key),
    enumerable: true,
  });
}

describe('webgpu device limits', () => {
  it('requests every maximum of the adapter', () => {
    const limits = getAdapterMaxLimits({ limits: new Limits() as unknown as GPUSupportedLimits });
    expect(limits).toEqual({
      maxTextureDimension2D: 16384,
      maxStorageBufferBindingSize: 2147483644,
      maxComputeWorkgroupsPerDimension: 65535,
    });
  });

  it('keeps the alignment limits at their defaults', () => {
    const limits = getAdapterMaxLimits({ limits: new Limits() as unknown as GPUSupportedLimits });
    expect(Object.keys(limits).some((k) => k.startsWith('min'))).toBe(false);
  });
});
