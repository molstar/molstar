/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  BufferTarget,
  canEncodeVideo,
  Mp4OutputFormat,
  WebMOutputFormat,
  MovOutputFormat,
  MkvOutputFormat,
  Output,
  VideoSample,
  VideoSampleSource,
  Quality,
  type VideoCodec,
  type VideoEncodingConfig,
  type QualityLevel,
} from 'mediabunny';
import { ParamDefinition as PD } from '@molstar/core/util/param-definition';
import type { VideoRendererBackend, VideoRendererConfig } from '../backend.js';

const formats = {
  mp4: () => new Mp4OutputFormat(),
  webm: () => new WebMOutputFormat(),
  mov: () => new MovOutputFormat(),
  mkv: () => new MkvOutputFormat(),
};
const labels: Record<VideoCodec, string> = {
  avc: 'H.264',
  hevc: 'H.265',
  vp8: 'VP8',
  vp9: 'VP9',
  av1: 'AV1',
  prores: 'ProRes',
};

function encodingConfig(values: Record<string, any>): VideoEncodingConfig {
  const { name: codec, params: p } = values.codec;
  const rate = p.rate;
  const quality =
    rate.name === 'bitrate'
      ? new Quality({ bitrate: Math.round(rate.params.bitrate * 1e6), bitrateMode: rate.params.mode })
      : new Quality({ quality: rate.name as QualityLevel, preferBitrate: true });
  return {
    codec,
    quality,
    keyFrameInterval: p.keyFrameInterval,
    hardwareAcceleration: p.hardwareAcceleration,
    fullCodecString: p.fullCodecString || undefined,
    alpha: p.alpha ?? 'discard',
    latencyMode: 'quality',
  };
}

function codecParams(codec: VideoCodec, container: string): PD.Params {
  const params: PD.Params = {
    rate: PD.MappedStatic(
      'high',
      {
        low: PD.EmptyGroup(),
        medium: PD.EmptyGroup(),
        high: PD.EmptyGroup(),
        'very-high': PD.EmptyGroup(),
        bitrate: PD.Group(
          {
            bitrate: PD.Numeric(
              8,
              { min: 0.1, max: 200, step: 0.1 },
              {
                label: 'Bitrate (Mbps)',
                help: () => ({
                  description:
                    'Target amount of video data per second, in megabits. Higher values generally preserve more detail and produce larger files. The actual bitrate depends on the codec and video content.',
                }),
              },
            ),
            mode: PD.Select(
              'variable',
              [
                ['variable', 'Variable'],
                ['constant', 'Constant'],
              ],
              {
                label: 'Bitrate mode',
                help: () => ({
                  description:
                    'Variable bitrate lets the encoder use more data for complex frames and less for simple frames. Constant bitrate aims to keep the data rate steady. Variable is the default for exported animations.',
                }),
              },
            ),
          },
          { isFlat: true },
        ),
      },
      {
        label: 'Quality',
        options: [
          ['low', 'Low'],
          ['medium', 'Medium'],
          ['high', 'High'],
          ['very-high', 'Very high'],
          ['bitrate', 'Custom bitrate'],
        ],
        help: () => ({
          description:
            'Quality presets let the encoder choose a bitrate for the codec, resolution and frame rate. Higher quality usually produces larger files. Choose Custom bitrate to set the target data rate yourself; High is the default.',
        }),
      },
    ),
    hardwareAcceleration: PD.Select(
      'no-preference',
      [
        ['no-preference', 'Automatic'],
        ['prefer-hardware', 'Prefer hardware'],
        ['prefer-software', 'Prefer software'],
      ],
      {
        label: 'Hardware acceleration',
        category: 'Advanced',
        help: () => ({
          description:
            'Automatic lets the encoder choose how to encode. Hardware encoding can be faster and use less CPU; software encoding uses the CPU. These are preferences, not guarantees, and availability depends on the selected codec and device.',
        }),
      },
    ),
    fullCodecString: PD.Text('', {
      label: 'Codec profile',
      placeholder: 'Automatic',
      category: 'Advanced',
      description: 'Leave empty for automatic selection.',
      help: () => ({
        description:
          'An optional full codec string specifying a profile and level for the selected codec. This can constrain compression features, resolution and playback compatibility. Leave empty to let the encoder choose; an unsupported or mismatched string prevents rendering.',
      }),
    }),
  };
  if (codec !== 'prores') {
    params.keyFrameInterval = PD.Numeric(
      2,
      { min: 0.1, max: 60, step: 0.1 },
      {
        label: 'Keyframe interval (s)',
        category: 'Advanced',
        help: () => ({
          description:
            'Time in seconds between independently encoded frames. Shorter intervals improve seeking but generally increase file size. Longer intervals allow more compression. The default is 2 seconds.',
        }),
      },
    );
  }
  if (codec === 'vp9' && (container === 'webm' || container === 'mkv')) {
    params.alpha = PD.Select(
      'discard',
      [
        ['discard', 'Discard'],
        ['keep', 'Keep'],
      ],
      {
        label: 'Transparency',
        category: 'Advanced',
        help: () => ({
          description:
            'Keep preserves transparency from the rendered frames when the encoder supports it. Discard exports opaque video. This option is available for VP9 in WebM or MKV and is only useful when the rendered frames contain transparent regions.',
        }),
      },
    );
  }
  return params;
}

/** Includes custom encoders registered with Mediabunny, including server-side ProRes encoders. */
export const MediabunnyBackend: VideoRendererBackend = {
  name: 'Mediabunny',
  dimensionAlignment: 1,
  async getParams(config, values = {}) {
    const containers = Object.entries(formats).map(([name, create]) => ({
      name,
      codecs: create().getSupportedVideoCodecs(),
    }));
    // Discover each codec once, then reuse its availability across containers.
    const codecs = [...new Set(containers.flatMap((format) => format.codecs))];
    const supported = new Set(
      (
        await Promise.all(
          codecs.map(async (codec) => {
            const params = codecParams(codec, values.container ?? 'mp4');
            const defaults = PD.getDefaultValues(params);
            const candidates = [defaults];
            if (values.codec?.name === codec) candidates.push(PD.normalizeParams(params, values.codec.params, 'all'));
            // A codec may support lower bitrates but not the High default (or require an explicit profile).
            for (const quality of ['low', 'medium', 'very-high']) {
              candidates.push({ ...defaults, rate: { name: quality, params: {} } });
            }
            for (const settings of candidates) {
              try {
                const encoding = encodingConfig({ codec: { name: codec, params: settings } });
                if (
                  await canEncodeVideo(codec, {
                    ...encoding,
                    width: config.width,
                    height: config.height,
                    frameRate: config.fps,
                  })
                )
                  return codec;
              } catch {
                // An invalid selected setting must not hide the schema needed to edit it; validation reports it separately.
              }
            }
            return undefined;
          }),
        )
      ).filter((codec) => codec !== undefined),
    );
    const availability = containers.map((format) => {
      const codecs = format.codecs.filter((codec) => supported.has(codec));
      // Preserve an existing compatible selection when FPS/resolution makes its encoder unavailable.
      // Validation disables rendering instead of silently replacing the codec or format.
      const current = values.codec?.name as VideoCodec | undefined;
      if (
        format.name === (values.container ?? 'mp4') &&
        current &&
        format.codecs.includes(current) &&
        !codecs.includes(current)
      ) {
        codecs.push(current);
      }
      return { ...format, codecs };
    });
    const available = availability.filter((f) => f.codecs.length);
    if (!available.length)
      throw new Error(
        'No video encoder is available. Use a browser with WebCodecs or explicitly provide a WASM/custom backend.',
      );
    const selected = available.find((f) => f.name === values.container) ?? available[0];
    const groups: PD.Params = {};
    for (const codec of selected.codecs) groups[codec] = PD.Group(codecParams(codec, selected.name), { isFlat: true });
    const defaultCodec = ['avc', 'vp8', 'vp9', 'hevc', 'av1', 'prores'].find((codec) =>
      selected.codecs.includes(codec as VideoCodec),
    )!;
    return {
      container: PD.Select(
        selected.name,
        available.map((f) => [f.name, f.name.toUpperCase()]),
        {
          label: 'Format',
          description: 'File format of the exported video.',
          help: () => ({
            description:
              'The format determines the file extension and which codecs can be used. The codec determines how the video is compressed. Changing the format updates the codec choices. Available encoders are checked at the current resolution and FPS. An unavailable current selection is retained so you can adjust its settings without changing the codec automatically.',
          }),
        },
      ),
      codec: PD.MappedStatic(defaultCodec, groups, {
        options: selected.codecs.map((codec) => [
          codec,
          `${labels[codec]}${supported.has(codec) ? '' : ' (unavailable)'}`,
        ]),
        help: () => ({
          description:
            'The codec compresses the video frames. Its settings control quality, file size and encoding performance. Choose a quality preset or a target bitrate; a higher bitrate generally preserves more detail and produces a larger file. Render is enabled once the selected settings have been checked for encoder support.',
        }),
      }),
    };
  },
  getPresets(params, values = {}) {
    const defaults = PD.getDefaultValues(params);
    const current = PD.normalizeParams(params, values, 'all');
    const preset = (quality: QualityLevel, fps: number) =>
      PD.normalizeParams(
        params,
        {
          ...defaults,
          container: current.container,
          fps,
          codec: { name: current.codec.name, params: { rate: { name: quality, params: {} } } },
        },
        'all',
      );
    return [
      [preset('high', 30), 'Balanced (default)'],
      [preset('medium', 15), 'Fast preview'],
      [preset('very-high', 30), 'High quality'],
      [preset('high', 60), 'Smooth motion'],
    ];
  },
  async validate(config: VideoRendererConfig, values) {
    const format = formats[values.container as keyof typeof formats]?.();
    if (!format || !format.getSupportedVideoCodecs().includes(values.codec?.name))
      return 'Unsupported codec/container combination.';
    const p = values.codec.params;
    if (p.rate.name === 'bitrate' && (!Number.isFinite(p.rate.params.bitrate) || p.rate.params.bitrate <= 0))
      return 'Bitrate must be positive.';
    if (p.keyFrameInterval !== undefined && (!Number.isFinite(p.keyFrameInterval) || p.keyFrameInterval <= 0))
      return 'Keyframe interval must be positive.';
    const encoding = encodingConfig(values);
    if (
      !(await canEncodeVideo(encoding.codec, {
        ...encoding,
        width: config.width,
        height: config.height,
        frameRate: config.fps,
      }))
    ) {
      return 'The selected codec/settings cannot be encoded at this resolution and frame rate.';
    }
  },
  async create(config, values) {
    const format = formats[values.container as keyof typeof formats]();
    const target = new BufferTarget();
    const output = new Output({ format, target });
    const source = new VideoSampleSource(encodingConfig(values));
    try {
      output.addVideoTrack(source, { frameRate: config.fps });
      await output.start();
    } catch (e) {
      await output.cancel();
      throw e;
    }
    return {
      async addFrame(data, timestamp, duration) {
        const sample = new VideoSample(data, {
          format: 'RGBA',
          codedWidth: config.width,
          codedHeight: config.height,
          timestamp,
          duration,
        });
        try {
          await source.add(sample);
        } finally {
          sample.close();
        }
      },
      async finalize() {
        await output.finalize();
        return {
          data: new Uint8Array(target.buffer!),
          mimeType: await output.getMimeType(),
          extension: values.container,
        };
      },
      async dispose() {
        if (output.state !== 'finalized' && output.state !== 'canceled') await output.cancel();
      },
    };
  },
};
