/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Node has no WebCodecs, so Mediabunny has no encoder here. Its own extension point (registerEncoder) takes a stand-in
// that "encodes" by keeping the samples (as 16-bit bytes in the packets), which is enough to run the real sample feed
// and the real MP4 and WebM muxers, and to check the output's container, type and size. The codecs are not tested: the
// browser's are. The stand-ins stay registered for this file only (each Jest test file has its own modules).

import { CustomAudioEncoder, EncodedPacket, type AudioSample, registerEncoder } from 'mediabunny';
import { mediabunnyEncoder, type MediabunnyCodec } from '../mediabunny.js';
import { text, tone } from './helpers.js';

const supports = new Set<string>();
const seen: { codec: string; sampleRate: number; bitrate: number | undefined; frames: number }[] = [];

class StandIn extends CustomAudioEncoder {
  private frames = 0;
  static override supports(codec: string): boolean {
    return supports.has(codec);
  }
  async init(): Promise<void> {}
  async encode(sample: AudioSample): Promise<void> {
    const bytes = new Uint8Array(sample.allocationSize({ planeIndex: 0, format: 's16' }));
    sample.copyTo(bytes, { planeIndex: 0, format: 's16' });
    this.frames += sample.numberOfFrames;
    // OpusHead / AudioSpecificConfig: what the containers need to describe the track.
    const description =
      this.codec === 'opus'
        ? new Uint8Array([0x4f, 0x70, 0x75, 0x73, 0x48, 0x65, 0x61, 0x64, 1, 1, 0, 0, 0x80, 0xbb, 0, 0, 0, 0, 0])
        : new Uint8Array([0x13, 0x08]);
    this.onPacket(new EncodedPacket(bytes, 'key', sample.timestamp, sample.duration), {
      decoderConfig: {
        codec: this.codec === 'opus' ? 'opus' : 'mp4a.40.2',
        numberOfChannels: 1,
        sampleRate: this.config.sampleRate,
        description,
      },
    });
  }
  async flush(): Promise<void> {
    seen.push({
      codec: this.codec,
      sampleRate: this.config.sampleRate,
      bitrate: this.config.bitrate,
      frames: this.frames,
    });
  }
  async close(): Promise<void> {}
}
registerEncoder(StandIn);

const wanted = (codecs: MediabunnyCodec[]) => {
  supports.clear();
  for (const c of codecs) supports.add(c);
  seen.length = 0;
};

describe('mediabunnyEncoder: the containers (stand-in encoders)', () => {
  it('AAC in an .m4a: audio/mp4, an MP4 file, every sample fed, at the bitrate', async () => {
    wanted(['aac', 'opus']);
    const pcm = tone(24000, 2.5);
    const audio = await mediabunnyEncoder().encode({ kind: 'pcm', pcm, sampleRate: 24000 });
    expect(audio.mimeType).toBe('audio/mp4');
    expect(audio.extension).toBe('m4a');
    expect(audio.durationMs).toBe(2500);
    expect(text(audio.bytes, 4)).toBe('ftyp');
    expect(audio.bytes.length).toBeGreaterThan(pcm.length * 2);
    expect(seen).toEqual([{ codec: 'aac', sampleRate: 24000, bitrate: 64_000, frames: pcm.length }]);
  });

  it('the stand-in is fed the samples as they are (the chunks do not shift them)', async () => {
    wanted(['aac']);
    const pcm = tone(8000, 1.2);
    const audio = await mediabunnyEncoder({ prefer: ['aac'] }).encode({ kind: 'pcm', pcm, sampleRate: 8000 });
    // mdat holds the packets in order: the samples, little-endian.
    const want = new Uint8Array(pcm.buffer, pcm.byteOffset, pcm.byteLength);
    const hay = Buffer.from(audio.bytes);
    expect(hay.indexOf(Buffer.from(want.subarray(0, 4000)))).toBeGreaterThan(0);
    expect(hay.indexOf(Buffer.from(want.subarray(want.length - 400)))).toBeGreaterThan(0);
  });

  it('Opus in WebM when AAC cannot be: audio/webm, a Matroska file', async () => {
    wanted(['opus']);
    const audio = await mediabunnyEncoder({ bitrate: 32_000 }).encode({
      kind: 'pcm',
      pcm: tone(24000, 1),
      sampleRate: 24000,
    });
    expect(audio.mimeType).toBe('audio/webm');
    expect(audio.extension).toBe('webm');
    expect(Array.from(audio.bytes.subarray(0, 4))).toEqual([0x1a, 0x45, 0xdf, 0xa3]);
    expect(seen).toMatchObject([{ codec: 'opus', bitrate: 32_000 }]);
  });

  it('the order is an option: Opus first when both can be', async () => {
    wanted(['aac', 'opus']);
    const audio = await mediabunnyEncoder({ prefer: ['opus', 'aac'] }).encode({
      kind: 'pcm',
      pcm: tone(24000, 0.5),
      sampleRate: 24000,
    });
    expect(audio.mimeType).toBe('audio/webm');
  });

  it('neither: WAV', async () => {
    wanted([]);
    const audio = await mediabunnyEncoder().encode({ kind: 'pcm', pcm: tone(24000, 0.5), sampleRate: 24000 });
    expect(audio.mimeType).toBe('audio/wav');
    expect(seen).toEqual([]);
  });

  it('an abort during the feed stops it, with nothing written', async () => {
    wanted(['aac']);
    const stop = new AbortController();
    const encoder = mediabunnyEncoder({
      canEncode: async () => {
        stop.abort();
        return true;
      },
    });
    await expect(
      encoder.encode({ kind: 'pcm', pcm: tone(24000, 3), sampleRate: 24000 }, stop.signal),
    ).rejects.toThrow();
    expect(seen).toEqual([]);
  });
});
