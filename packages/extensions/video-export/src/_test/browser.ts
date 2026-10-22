/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { MediabunnyBackend } from '../backends/mediabunny.js';
import { encodeVideoAnimation } from '../encoder.js';
import { VideoExport } from '../index.js';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import { RuntimeContext } from '@molstar/core/task';
import { Input, BufferSource, ALL_FORMATS, EncodedPacketSink } from 'mediabunny';

export { VideoExport };

export async function run() {
  const config = { width: 64, height: 64, fps: 10 };
  const schema = await MediabunnyBackend.getParams(config);
  const results: string[] = [];
  const noop = async () => {};
  const plugin = {
    animationLoop: { isAnimating: false, stop() {}, resetTime() {}, tick: noop },
    managers: { animation: { play: noop, stop: noop } },
  };
  for (const [container] of (schema.container as PD.Select<string>).options) {
    const params = await MediabunnyBackend.getParams(config, { container });
    for (const [codec] of (params.codec as PD.Mapped<any>).select.options) {
      const values = PD.normalizeParams(params, { container, codec: { name: codec, params: {} } }, 'all');
      if (container === 'webm' && !['vp8', 'vp9', 'av1'].includes(codec)) throw new Error('Invalid WebM codec');
      if (codec === 'vp9' && container === 'webm' && !values.codec.params.alpha)
        throw new Error('Missing transparency param');
      if (container === 'mp4' && values.codec.params.alpha) throw new Error('Unexpected transparency param');
      if (results.length % 2) {
        values.codec.params.rate = { name: 'bitrate', params: { bitrate: 8.3, mode: 'variable' } };
      }
      const video = await encodeVideoAnimation(plugin as any, RuntimeContext.Synchronous, {
        ...config,
        viewport: { x: 0, y: 0, width: 64, height: 64 },
        animation: { definition: {} as any, params: {}, customDurationMs: 250 },
        backendParams: values,
        pass: { updateBackground: noop, getImageData: async () => ({ data: new Uint8Array(64 * 64 * 4).fill(255) }) },
      });
      if (video.extension !== container || !video.mimeType.startsWith('video/'))
        throw new Error('Invalid video metadata');
      const input = new Input({ source: new BufferSource(video.data), formats: ALL_FORMATS });
      try {
        const track = await input.getPrimaryVideoTrack();
        if (!track || track.codec !== codec) throw new Error('Incorrect encoded codec');
        let count = 0;
        for await (const _packet of new EncodedPacketSink(track).packets()) count++;
        if (count !== 3) throw new Error(`Expected 3 frames, got ${count}`);
        const duration = await input.computeDuration();
        if (duration < 0.24 || duration > 0.31) throw new Error(`Invalid duration: ${duration}`);
      } finally {
        input.dispose();
      }
      results.push(`${container}/${codec}`);
    }
  }
  // Selecting a preset must preserve a non-default codec and its format.
  const presetSchema = { fps: PD.Numeric(30), ...(await MediabunnyBackend.getParams(config, { container: 'webm' })) };
  const codec = (presetSchema.codec as PD.Mapped<any>).select.options.at(-1)![0];
  const current = PD.normalizeParams(presetSchema, { container: 'webm', codec: { name: codec, params: {} } }, 'all');
  const presets = MediabunnyBackend.getPresets!(presetSchema, current);
  for (const [preset, name] of presets) {
    if (preset.container !== 'webm' || preset.codec.name !== codec)
      throw new Error(`Preset changed format/codec: ${name}`);
    if (preset.codec.params.hardwareAcceleration !== 'no-preference')
      throw new Error(`Preset changed Automatic: ${name}`);
    const video = await encodeVideoAnimation(plugin as any, RuntimeContext.Synchronous, {
      ...config,
      fps: preset.fps,
      viewport: { x: 0, y: 0, width: 64, height: 64 },
      animation: { definition: {} as any, params: {}, customDurationMs: 250 },
      backendParams: preset,
      pass: { updateBackground: noop, getImageData: async () => ({ data: new Uint8Array(64 * 64 * 4).fill(255) }) },
    });
    if (!video.data.length) throw new Error(`Preset produced no video: ${name}`);
    results.push(`preset/${name}`);
  }
  return results;
}
