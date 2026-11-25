/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Tadej Satler <tadej.satler@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { ValueCell } from '@molstar/core/util';
import { DirectVolumeRenderableSchema } from '@molstar/graphics/gl/renderable/direct-volume';
import { DefineSpec } from '@molstar/graphics/gl/renderable/schema';
import { createWebGPUContext, requestWebGPUDevice } from '@molstar/graphics/gl/webgpu/context';
import { getDataLayerSchema } from '@molstar/graphics/gl/webgpu/data-layers';
import { createWebGPUProgram } from '@molstar/graphics/gl/webgpu/program';
import { getWgslProgram, resolveWgslSource } from '@molstar/graphics/gl/webgpu/shader/programs';
import { WebGPU } from './webgpu-module.js';

describe.skipIf(!WebGPU.module)('direct-volume WGSL', () => {
  it('compiles in every render mode and variant', async () => {
    const gpu = WebGPU.module!.create([]);
    // a copy: the data layer schema is shared (cached per schema)
    const schema = { ...getDataLayerSchema(DirectVolumeRenderableSchema) };
    Object.defineProperty(schema, 'dRenderVariant', { value: DefineSpec('string') });
    const base: Record<string, unknown> = {};
    for (const [k, s] of Object.entries(schema)) {
      if (s.type === 'define') base[k] = s.kind === 'string' ? s.options?.[0] : s.kind === 'boolean' ? false : 1;
    }
    const errors: string[] = [];
    for (const dRenderMode of ['composite', 'maximum', 'attenuated', 'additive']) {
      for (const dRenderVariant of ['colorBlended', 'colorWboit', 'colorDpoit', 'marking', 'pick']) {
        for (const dClipObjects of [false, true]) {
          // a device per program: Dawn's Node bindings crash compiling several of these on one device
          const device = await requestWebGPUDevice('default', gpu);
          if (!device) throw new Error('no WebGPU device');
          const graphics = createWebGPUContext(device, { width: 16, height: 16, ownsDevice: true });
          try {
            const v = {
              ...base,
              dGridTexType: '3d',
              dColorType: 'uniform',
              dLightCount: 2,
              dRenderMode,
              dRenderVariant,
              dClipObjects,
              dClipVariant: 'pixel',
            };
            const cells = Object.fromEntries(Object.entries(v).map(([k, x]) => [k, ValueCell.create(x)]));
            const source = resolveWgslSource(getWgslProgram('direct-volume')!.source, cells);
            const p = createWebGPUProgram(device, 'direct-volume', source, cells, schema, (m) =>
              errors.push(`${dRenderMode} ${dRenderVariant}${dClipObjects ? ' clipped' : ''}: ${m}`),
            );
            await p.module.getCompilationInfo();
          } finally {
            graphics.destroy();
          }
        }
      }
    }
    expect(errors).toEqual([]);
  });
});
