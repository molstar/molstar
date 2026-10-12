/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { BehaviorSubject, Subject } from 'rxjs';
import { VideoControls } from '../controls.js';
import { CustomVideoEncoder, registerEncoder, type VideoCodec } from 'mediabunny';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { resolveVideoParams } from '../backend.js';
import { MediabunnyBackend } from '../backends/mediabunny.js';

// Use real Mediabunny capability checks without depending on browser/platform codecs.
// Each mode uses different dimensions because Mediabunny memoizes capability results.
class TestEncoder extends CustomVideoEncoder {
  static mode: 'all' | 'low' | 'profile' | 'fps' = 'all';
  static calls: { codec: VideoCodec; config: VideoEncoderConfig }[] = [];
  static supports(codec: VideoCodec, config: VideoEncoderConfig) {
    this.calls.push({ codec, config });
    if (this.mode === 'low') return codec === 'vp8' && config.bitrate! < 20000;
    if (this.mode === 'fps') return codec === 'avc' || (codec === 'av1' && config.framerate! <= 30);
    if (this.mode === 'profile') return codec === 'hevc' && config.codec === 'hvc1.2.4.L123.B0';
    return true;
  }
  init() {}
  encode() {}
  flush() {}
  close() {}
}
registerEncoder(TestEncoder);
beforeEach(() => {
  TestEncoder.mode = 'all';
  TestEncoder.calls = [];
});

function values(codec = 'avc', bitrate?: number) {
  return {
    container: 'mp4',
    codec: {
      name: codec,
      params: {
        rate:
          bitrate === undefined
            ? { name: 'high', params: {} }
            : { name: 'bitrate', params: { bitrate, mode: 'variable' } },
        hardwareAcceleration: 'no-preference',
        keyFrameInterval: 2,
        fullCodecString: '',
      },
    },
  };
}

it.each([8.3, 16.6])('accepts fractional Mbps and converts it to integer bits per second', async (bitrate) => {
  const input = values('avc', bitrate);
  expect(await MediabunnyBackend.validate({ width: 64, height: 64, fps: 30 }, input)).toBeUndefined();
  expect(TestEncoder.calls.at(-1)!.config.bitrate).toBe(Math.round(bitrate * 1e6));
});

it.each([0, -1, NaN, Infinity])('reports invalid bitrates before constructing Quality', async (bitrate) => {
  expect(await MediabunnyBackend.validate({ width: 64, height: 64, fps: 30 }, values('avc', bitrate))).toMatch(
    /Bitrate/,
  );
  expect(TestEncoder.calls).toHaveLength(0);
});

it('reuses one codec discovery across all compatible formats', async () => {
  const params = await MediabunnyBackend.getParams({ width: 128, height: 128, fps: 30 });
  expect((params.container as PD.Select<string>).options.map(([name]) => name)).toEqual(['mp4', 'webm', 'mov', 'mkv']);
  expect(TestEncoder.calls).toHaveLength(6);
  expect(new Set(TestEncoder.calls.map((call) => call.codec)).size).toBe(6);
  const defaults = PD.getDefaultValues(params);
  expect(defaults.codec.params.hardwareAcceleration).toBe('no-preference');
});

it('discovers and resolves an encoder which supports Low but not High quality', async () => {
  TestEncoder.mode = 'low';
  const config = { width: 96, height: 96, fps: 30 };
  const params = await MediabunnyBackend.getParams(config);
  expect((params.codec as PD.Mapped<any>).select.options.map(([codec]) => codec)).toEqual(['vp8']);
  const input = { container: 'webm', codec: { name: 'vp8', params: { rate: { name: 'low', params: {} } } } };
  const normalized = await resolveVideoParams(MediabunnyBackend, config, input);
  expect(normalized.codec.params.rate.name).toBe('low');
});

it('includes a requested codec supported only with its explicit profile', async () => {
  TestEncoder.mode = 'profile';
  const input = values('hevc');
  input.codec.params.fullCodecString = 'hvc1.2.4.L123.B0';
  const normalized = await resolveVideoParams(MediabunnyBackend, { width: 160, height: 160, fps: 30 }, input);
  expect(normalized.codec.params.fullCodecString).toBe('hvc1.2.4.L123.B0');
});

it.each([
  ['mp4', 'av1'],
  ['webm', 'vp9'],
])('presets preserve the selected format/codec and retain Automatic', async (container, codec) => {
  const params = {
    fps: PD.Numeric(30),
    ...(await MediabunnyBackend.getParams({ width: 192, height: 192, fps: 30 }, { container })),
  };
  const current = PD.normalizeParams(
    params,
    { container, codec: { name: codec, params: { hardwareAcceleration: 'prefer-software' } } },
    'all',
  );
  const presets = MediabunnyBackend.getPresets!(params, current);
  expect(presets.map(([value]) => [value.fps, value.codec.params.rate.name])).toEqual([
    [30, 'high'],
    [15, 'medium'],
    [30, 'very-high'],
    [60, 'high'],
  ]);
  for (const [preset] of presets) {
    expect(preset.container).toBe(container);
    expect(preset.codec.name).toBe(codec);
    expect(preset.codec.params.hardwareAcceleration).toBe('no-preference');
  }
});

it.each(['mp4', 'webm'])('keeps the selected format/codec when a preset FPS is unsupported', async (container) => {
  TestEncoder.mode = 'fps';
  const plugin = {
    canvas3d: { resized: new Subject<void>() },
    helpers: {
      viewportScreenshot: {
        getSizeAndViewport: () => ({ viewport: { width: 224, height: 224 } }),
        events: { previewed: new Subject<void>() },
      },
    },
    managers: {
      animation: { animations: [], events: { updated: new Subject<void>() } },
      snapshot: { events: { changed: new Subject<void>() } },
    },
    behaviors: { state: { isBusy: new BehaviorSubject(false) } },
  };
  const controls = new VideoControls(plugin as any, MediabunnyBackend);
  try {
    await new Promise(setImmediate);
    controls.setParams({ ...controls.behaviors.params.value, container, codec: { name: 'av1', params: {} } });
    await new Promise(setImmediate);
    expect(controls.behaviors.encoder.value.canApply).toBe(true);
    const smooth = controls.videoParams.presets!.find(([, name]) => name === 'Smooth motion')![0];
    controls.setParams(smooth);
    await new Promise(setImmediate);
    expect(controls.behaviors.params.value.codec.name).toBe('av1');
    expect(controls.behaviors.params.value.container).toBe(container);
    expect(controls.behaviors.params.value.fps).toBe(60);
    expect(controls.behaviors.encoder.value.canApply).toBe(false);
    controls.setParams({ ...controls.behaviors.params.value, fps: 30 });
    await new Promise(setImmediate);
    expect(controls.behaviors.params.value.codec.name).toBe('av1');
    expect(controls.behaviors.params.value.container).toBe(container);
    expect(controls.behaviors.encoder.value.canApply).toBe(true);
  } finally {
    controls.dispose();
  }
});
