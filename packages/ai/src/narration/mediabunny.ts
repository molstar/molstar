/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The narration encoder over Mediabunny (.v6/ai/plan.md): the browser's own AAC and Opus encoders (WebCodecs), found by
 * asking whether it can encode, never by which browser it is. At its own entry (`./mediabunny`) so the core entry does
 * not pull Mediabunny in. No MP3 is encoded: its encoders are LGPL, and a provider's MP3 is kept as it arrives.
 */

import {
  AudioSample,
  AudioSampleSource,
  BufferTarget,
  canEncodeAudio,
  Mp4OutputFormat,
  Output,
  type OutputFormat,
  Quality,
  WebMOutputFormat,
} from 'mediabunny';
import { pcmDuration } from './audio.js';
import { type NarrationAudio, type NarrationEncoder, passThrough, wavEncoder } from './encoder.js';

export type MediabunnyCodec = 'aac' | 'opus';

export interface MediabunnyEncoderOptions {
  /** The codecs to try, first choice first (the first the browser can encode is used). Default: AAC, then Opus. */
  prefer?: readonly MediabunnyCodec[];
  /** Bits per second. Default 64 000: speech. */
  bitrate?: number;
  /** Whether a codec can be encoded here, for this audio. Default: Mediabunny's check (WebCodecs'). */
  canEncode?: (codec: MediabunnyCodec, audio: { sampleRate: number; bitrate: number }) => Promise<boolean>;
}

const DEFAULT_BITRATE = 64_000;

/** About a second of samples at a time: a WebCodecs encoder takes a stream, not a minute of sound at once. */
const CHUNK_SECONDS = 1;

const FORMATS: Record<MediabunnyCodec, { format: () => OutputFormat; mimeType: string; extension: string }> = {
  // The metadata goes first in the .m4a, so it plays as it loads.
  aac: { format: () => new Mp4OutputFormat({ fastStart: 'in-memory' }), mimeType: 'audio/mp4', extension: 'm4a' },
  opus: { format: () => new WebMOutputFormat(), mimeType: 'audio/webm', extension: 'webm' },
};

/**
 * PCM as AAC in an .m4a (audio/mp4) when the browser can encode AAC, else Opus in WebM, else WAV (the order is the
 * `prefer` option). Audio a provider already encoded passes through unchanged.
 */
export function mediabunnyEncoder(options: MediabunnyEncoderOptions = {}): NarrationEncoder {
  const prefer = options.prefer ?? ['aac', 'opus'];
  const bitrate = options.bitrate ?? DEFAULT_BITRATE;
  const canEncode =
    options.canEncode ??
    ((codec, { sampleRate }) =>
      canEncodeAudio(codec, { numberOfChannels: 1, sampleRate, quality: new Quality({ bitrate }) }));
  const wav = wavEncoder();
  return {
    async encode(audio, signal) {
      signal?.throwIfAborted();
      if (audio.kind === 'encoded') return passThrough(audio);
      for (const codec of prefer) {
        if (!(await supported(canEncode, codec, audio.sampleRate, bitrate))) continue;
        signal?.throwIfAborted();
        try {
          return await encodeWith(codec, audio.pcm, audio.sampleRate, bitrate, signal);
        } catch (error) {
          // An encoder that said it could and did not (or a container that would not form): the next one, then WAV.
          if (signal?.aborted) throw error;
        }
      }
      return wav.encode(audio, signal);
    },
  };
}

async function supported(
  canEncode: NonNullable<MediabunnyEncoderOptions['canEncode']>,
  codec: MediabunnyCodec,
  sampleRate: number,
  bitrate: number,
): Promise<boolean> {
  try {
    return await canEncode(codec, { sampleRate, bitrate });
  } catch {
    // The check itself failed (an unusual rate, no WebCodecs): it cannot.
    return false;
  }
}

async function encodeWith(
  codec: MediabunnyCodec,
  pcm: Int16Array,
  sampleRate: number,
  bitrate: number,
  signal: AbortSignal | undefined,
): Promise<NarrationAudio> {
  const { format, mimeType, extension } = FORMATS[codec];
  const target = new BufferTarget();
  const output = new Output({ format: format(), target });
  const source = new AudioSampleSource({ codec, quality: new Quality({ bitrate }) });
  output.addAudioTrack(source);
  try {
    await output.start();
    const chunk = Math.max(1, Math.round(sampleRate * CHUNK_SECONDS));
    for (let at = 0; at < pcm.length; at += chunk) {
      signal?.throwIfAborted();
      const part = pcm.subarray(at, at + chunk);
      const sample = new AudioSample({
        data: part,
        format: 's16',
        numberOfChannels: 1,
        sampleRate,
        timestamp: at / sampleRate,
      });
      try {
        await source.add(sample);
      } finally {
        sample.close();
      }
    }
    await output.finalize();
  } catch (error) {
    await output.cancel().catch(() => undefined);
    throw error;
  }
  const buffer = target.buffer;
  if (!buffer || buffer.byteLength === 0) throw new Error('The encoder gave no audio.');
  return { bytes: new Uint8Array(buffer), mimeType, extension, durationMs: pcmDuration(pcm, sampleRate) };
}
