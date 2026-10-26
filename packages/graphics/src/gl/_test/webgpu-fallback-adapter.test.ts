/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { describe, expect, it, vi } from 'vitest';
import { AutoWebGPUBackend } from '../webgpu/backend.js';
import { requestWebGPUDevice } from '../webgpu/context.js';

const device = { label: 'stub device' };

/** A `GPU` of one adapter with `float32-blendable`, a fallback (software) adapter or not */
function stubGpu(isFallbackAdapter: boolean) {
  const adapter = {
    info: { isFallbackAdapter },
    features: new Set(['float32-blendable']),
    limits: {},
    requestDevice: vi.fn(async () => device),
  };
  return { gpu: { requestAdapter: async () => adapter } as unknown as GPU, adapter };
}

describe('WebGPU fallback adapters', () => {
  // auto took SwiftShader's WebGPU over WebGL2 on the GPU: a software adapter has `float32-blendable` too
  it('are rejected where asked (the auto backend), used otherwise', async () => {
    const { gpu, adapter } = stubGpu(true);
    expect(await requestWebGPUDevice('default', gpu, 'adapter', 'none', [], true)).toBeUndefined();
    expect(adapter.requestDevice).not.toHaveBeenCalled();
    expect(await requestWebGPUDevice('default', gpu, 'adapter', 'none', [], false)).toBe(device);
  });

  it('create no context of the auto backend, so that the next backend of the list is used', async () => {
    const { gpu, adapter } = stubGpu(true);
    vi.stubGlobal('navigator', { gpu });
    expect(AutoWebGPUBackend.isAvailable()).toBe(true);
    expect(await AutoWebGPUBackend.create({ width: 1, height: 1 } as HTMLCanvasElement, {} as any)).toBeUndefined();
    expect(adapter.requestDevice).not.toHaveBeenCalled();
  });

  it('do not reject hardware adapters', async () => {
    const { gpu } = stubGpu(false);
    expect(await requestWebGPUDevice('default', gpu, 'adapter', 'none', [], true)).toBe(device);
  });
});
