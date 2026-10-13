/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { RuntimeContext } from '@molstar/core/task';
import { encodeVideoAnimation } from '../encoder.js';
import { resolveVideoParams, type VideoRendererBackend } from '../backend.js';
import { H264WasmBackend } from '../backends/h264-wasm.js';

function fixture() {
  const frames: { timestamp: number; duration: number }[] = [];
  const dispose = jest.fn();
  const backend: VideoRendererBackend = {
    name: 'Test',
    dimensionAlignment: 1,
    getParams: async () => ({}),
    validate: async () => undefined,
    create: async () => ({
      addFrame: (_data, timestamp, duration) => {
        frames.push({ timestamp, duration });
      },
      finalize: () => ({ data: new Uint8Array([1]), mimeType: 'video/test', extension: 'test' }),
      dispose,
    }),
  };
  const plugin = {
    animationLoop: {
      isAnimating: true,
      stop: jest.fn(),
      start: jest.fn(),
      resetTime: jest.fn(),
      tick: jest.fn(async () => {}),
    },
    managers: { animation: { play: jest.fn(async () => {}), stop: jest.fn(async () => {}) } },
    canvas3d: { props: { renderer: { backgroundColor: 0 }, transparentBackground: true }, setProps: jest.fn() },
  };
  const params = {
    animation: { definition: {} as any, params: {}, customDurationMs: 1250 },
    width: 16,
    height: 16,
    viewport: { x: 0, y: 0, width: 16, height: 16 },
    fps: 2,
    pass: {
      updateBackground: async () => {},
      getImageData: async () => ({ data: new Uint8Array(16 * 16 * 4) }),
    } as any,
    backend,
  };
  return { plugin, params, frames, dispose };
}

it('timestamps frames at the requested FPS and bounds the final sample duration', async () => {
  const f = fixture();
  const result = await encodeVideoAnimation(f.plugin as any, RuntimeContext.Synchronous, f.params);
  expect(result.extension).toBe('test');
  expect(f.frames).toEqual([
    { timestamp: 0, duration: 0.5 },
    { timestamp: 0.5, duration: 0.5 },
    { timestamp: 1, duration: 0.25 },
  ]);
  expect(f.plugin.animationLoop.tick.mock.calls).toHaveLength(3);
  expect(f.dispose).toHaveBeenCalledTimes(1);
  expect(f.plugin.animationLoop.start).toHaveBeenCalledTimes(1);
});

it('releases the backend and restores playback and background on failure/abort', async () => {
  const f = fixture();
  const error = new Error('Aborted');
  f.params.pass.getImageData = async () => {
    throw error;
  };
  await expect(
    encodeVideoAnimation(f.plugin as any, RuntimeContext.Synchronous, { ...f.params, customBackground: 1 as any }),
  ).rejects.toBe(error);
  expect(f.dispose).toHaveBeenCalledTimes(1);
  expect(f.plugin.managers.animation.stop).toHaveBeenCalledTimes(1);
  expect(f.plugin.canvas3d.setProps).toHaveBeenLastCalledWith({
    renderer: { backgroundColor: 0 },
    transparentBackground: true,
  });
  expect(f.plugin.animationLoop.start).toHaveBeenCalledTimes(1);
});

it('rejects invalid FPS before creating a backend', async () => {
  const f = fixture();
  f.params.backend.create = jest.fn(f.params.backend.create);
  await expect(
    encodeVideoAnimation(f.plugin as any, RuntimeContext.Synchronous, { ...f.params, fps: 0 }),
  ).rejects.toThrow('Frame rate');
  expect(f.params.backend.create).not.toHaveBeenCalled();
});

it('does not substitute an explicitly requested unavailable codec', async () => {
  const f = fixture();
  f.params.backend.getParams = async () => ({ codec: PD.MappedStatic('avc', { avc: PD.Group({}) }) });
  await expect(
    resolveVideoParams(
      f.params.backend,
      { width: 16, height: 16, fps: 30 },
      {
        codec: { name: 'av1', params: {} },
      },
    ),
  ).rejects.toThrow('Unsupported video codec');
});

it('validates WASM quantization and dimension requirements', async () => {
  expect(await H264WasmBackend.validate({ width: 16, height: 16, fps: 30 }, { quantization: 18 })).toBeUndefined();
  expect(await H264WasmBackend.validate({ width: 15, height: 16, fps: 30 }, { quantization: 18 })).toMatch(/even/);
  expect(await H264WasmBackend.validate({ width: 16, height: 16, fps: 30 }, { quantization: 52 })).toMatch(
    /Quantization/,
  );
});

it('provides WASM presets through PD with defaults matching Balanced', async () => {
  const params = { fps: PD.Numeric(30), ...(await H264WasmBackend.getParams({ width: 16, height: 16, fps: 30 })) };
  const presets = H264WasmBackend.getPresets!(params);
  const group = PD.Group(params, { presets });
  expect(group.presets![0][0]).toEqual(group.defaultValue);
  for (const [value] of presets) {
    expect(await H264WasmBackend.validate({ width: 16, height: 16, fps: value.fps }, value)).toBeUndefined();
  }
});
