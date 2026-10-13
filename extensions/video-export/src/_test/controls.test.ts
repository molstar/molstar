/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { expect, it, vi } from 'vitest';
import { BehaviorSubject, Subject } from 'rxjs';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { VideoControls } from '../controls.js';
import type { VideoRendererBackend } from '../backend.js';

function deferred<T>() {
  let resolve!: (value: T) => void;
  const promise = new Promise<T>((r) => {
    resolve = r;
  });
  return { promise, resolve };
}
const flush = async () => {
  await Promise.resolve();
  await Promise.resolve();
  await Promise.resolve();
};

function fixture(getParams: VideoRendererBackend['getParams'], overrides: Partial<VideoRendererBackend> = {}) {
  const backend: VideoRendererBackend = {
    name: 'Test',
    dimensionAlignment: 1,
    getParams,
    validate: async () => undefined,
    create: async () => {
      throw new Error('Not used');
    },
    ...overrides,
  };
  const resized = new Subject<void>();
  const plugin = {
    canvas3d: { resized },
    helpers: {
      viewportScreenshot: {
        getSizeAndViewport: () => ({ viewport: { width: 16, height: 16 } }),
        events: { previewed: new Subject<void>() },
      },
    },
    managers: {
      animation: { animations: [], events: { updated: new Subject<void>() } },
      snapshot: { events: { changed: new Subject<void>() } },
    },
    behaviors: { state: { isBusy: new BehaviorSubject(false) } },
  };
  return new VideoControls(plugin as any, backend);
}

it('ignores stale asynchronous schema results', async () => {
  const old = deferred<PD.Params>();
  const current = deferred<PD.Params>();
  const getParams = vi.fn().mockReturnValueOnce(old.promise).mockReturnValueOnce(current.promise);
  const controls = fixture(getParams);
  controls.setParams({ fps: 60, bitrate: 12 });
  expect(controls.behaviors.encoder.value.canApply).toBe(false);
  current.resolve({ bitrate: PD.Numeric(8) });
  await flush();
  expect(controls.behaviors.params.value).toEqual({ fps: 60, bitrate: 12 });
  expect(controls.behaviors.encoder.value.canApply).toBe(true);
  old.resolve({ quantization: PD.Numeric(18) });
  await flush();
  expect(controls.behaviors.paramsSchema.value.quantization).toBeUndefined();
  expect(controls.behaviors.params.value).toEqual({ fps: 60, bitrate: 12 });
  controls.dispose();
});

it('preserves valid values and replaces unavailable codec settings with new defaults', async () => {
  let codec = 'avc';
  const controls = fixture(async () => ({
    codec: PD.MappedStatic(codec, {
      [codec]: PD.Group(codec === 'avc' ? { quantization: PD.Numeric(18) } : { bitrate: PD.Numeric(8) }),
    }),
  }));
  await flush();
  controls.setParams({ fps: 60, codec: { name: 'avc', params: { quantization: 20 } } });
  await flush();
  expect(controls.behaviors.params.value.codec.params.quantization).toBe(20);
  codec = 'vp9';
  controls.setParams({ ...controls.behaviors.params.value, fps: 30 });
  await flush();
  expect(controls.behaviors.params.value).toEqual({ fps: 30, codec: { name: 'vp9', params: { bitrate: 8 } } });
  controls.dispose();
});

it('does not publish pending results after disposal', async () => {
  const pending = deferred<PD.Params>();
  const controls = fixture(() => pending.promise);
  const next = vi.spyOn(controls.behaviors.paramsSchema, 'next');
  controls.dispose();
  pending.resolve({ bitrate: PD.Numeric(8) });
  await flush();
  expect(next).not.toHaveBeenCalled();
});

it('validates ordinary edits without rediscovery, but regenerates for FPS, format and resolution', async () => {
  const getParams = vi.fn(async () => ({
    container: PD.Select('mp4', [
      ['mp4', 'MP4'],
      ['webm', 'WebM'],
    ]),
    quality: PD.Select('high', [
      ['high', 'High'],
      ['low', 'Low'],
    ]),
  }));
  const validate = vi.fn(async () => undefined);
  const controls = fixture(getParams, { validate });
  await flush();
  controls.setParams({ ...controls.behaviors.params.value, quality: 'low' });
  await flush();
  expect(getParams).toHaveBeenCalledTimes(1);
  expect(validate).toHaveBeenCalledTimes(2);
  controls.setParams({ ...controls.behaviors.params.value, fps: 60 });
  await flush();
  controls.setParams({ ...controls.behaviors.params.value, container: 'webm' });
  await flush();
  expect(getParams).toHaveBeenCalledTimes(3);
  const plugin = (controls as any).plugin;
  plugin.helpers.viewportScreenshot.getSizeAndViewport = () => ({ viewport: { width: 32, height: 32 } });
  plugin.canvas3d.resized.next();
  await flush();
  expect(getParams).toHaveBeenCalledTimes(4);
  expect(controls.behaviors.params.value.quality).toBe('low');
  controls.dispose();
});

it('ignores stale validation results when the schema is reused', async () => {
  const old = deferred<string | undefined>();
  const current = deferred<string | undefined>();
  const validate = vi
    .fn()
    .mockResolvedValueOnce(undefined)
    .mockReturnValueOnce(old.promise)
    .mockReturnValueOnce(current.promise);
  const getParams = vi.fn(async () => ({ bitrate: PD.Numeric(8) }));
  const controls = fixture(getParams, { validate });
  await flush();
  controls.setParams({ fps: 30, bitrate: 9 });
  controls.setParams({ fps: 30, bitrate: 10 });
  current.resolve(undefined);
  await flush();
  old.resolve('Unsupported old settings');
  await flush();
  expect(getParams).toHaveBeenCalledTimes(1);
  expect(controls.behaviors.encoder.value.canApply).toBe(true);
  expect(controls.behaviors.params.value.bitrate).toBe(10);
  controls.dispose();
});

it('exposes PD presets for a backend without codec/quantization parameter names', async () => {
  const presets = [[{ fps: 15, custom: 2 }, 'Preview']] as [Record<string, any>, string][];
  const getPresets = vi.fn(() => presets);
  const controls = fixture(async () => ({ custom: PD.Numeric(1) }), { getPresets });
  await flush();
  expect(controls.videoParams.presets).toBe(presets);
  expect(getPresets).toHaveBeenCalledWith(controls.behaviors.paramsSchema.value, controls.behaviors.params.value);
  controls.setParams(presets[0][0]);
  await flush();
  expect(controls.behaviors.params.value).toEqual({ fps: 15, custom: 2 });
  controls.dispose();
});
