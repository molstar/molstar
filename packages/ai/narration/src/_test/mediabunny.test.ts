/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import { type NarrationAudio, type ProviderAudio, wavSamples } from '../index.js';
import { type MediabunnyCodec, mediabunnyEncoder } from '../mediabunny.js';
import { mp3Frames, text, tone } from './helpers.js';

const pcm = (seconds = 0.5, sampleRate = 24000): ProviderAudio => ({
  kind: 'pcm',
  pcm: tone(sampleRate, seconds),
  sampleRate,
});

/** The codecs asked about, and an answer for each. */
function probe(can: Partial<Record<MediabunnyCodec, boolean>>) {
  const asked: [MediabunnyCodec, number, number][] = [];
  return {
    asked,
    canEncode: async (codec: MediabunnyCodec, { sampleRate, bitrate }: { sampleRate: number; bitrate: number }) => {
      asked.push([codec, sampleRate, bitrate]);
      return can[codec] ?? false;
    },
  };
}

const isWav = (a: NarrationAudio) => a.mimeType === 'audio/wav' && text(a.bytes, 0) === 'RIFF';

describe('mediabunnyEncoder: the choice', () => {
  it('in Node (no WebCodecs) it is a WAV', async () => {
    expect(typeof (globalThis as { AudioEncoder?: unknown }).AudioEncoder).toBe('undefined');
    const audio = await mediabunnyEncoder().encode(pcm());
    expect(isWav(audio)).toBe(true);
    expect(audio.durationMs).toBe(500);
    expect(wavSamples(audio.bytes)?.sampleRate).toBe(24000);
  });

  it("neither codec encodable: WAV, after asking for each in order, at the audio's rate and the bitrate", async () => {
    const p = probe({});
    const audio = await mediabunnyEncoder({ canEncode: p.canEncode, bitrate: 48_000 }).encode(pcm(0.2, 16000));
    expect(isWav(audio)).toBe(true);
    expect(p.asked).toEqual([
      ['aac', 16000, 48_000],
      ['opus', 16000, 48_000],
    ]);
  });

  it('a default bitrate suited to speech: 64 kbps', async () => {
    const p = probe({});
    await mediabunnyEncoder({ canEncode: p.canEncode }).encode(pcm());
    expect(p.asked[0]?.[2]).toBe(64_000);
  });

  it('the check that fails is a codec that cannot', async () => {
    const audio = await mediabunnyEncoder({
      canEncode: async () => {
        throw new TypeError('unusual rate');
      },
    }).encode(pcm());
    expect(isWav(audio)).toBe(true);
  });

  it('an MP3 passes through without asking', async () => {
    const p = probe({ aac: true });
    const bytes = mp3Frames(1);
    const audio = await mediabunnyEncoder({ canEncode: p.canEncode }).encode({
      kind: 'encoded',
      bytes,
      mimeType: 'audio/mpeg',
      durationMs: 1000,
    });
    expect(audio.bytes).toBe(bytes);
    expect(audio.extension).toBe('mp3');
    expect(p.asked).toEqual([]);
  });

  it('an abort before it starts stops it', async () => {
    const stop = new AbortController();
    stop.abort();
    await expect(mediabunnyEncoder().encode(pcm(), stop.signal)).rejects.toThrow();
  });

  it('a codec that said it could and could not falls to the next, then to WAV', async () => {
    // Without WebCodecs, Mediabunny cannot encode either, whatever the probe says.
    const p = probe({ aac: true, opus: true });
    const audio = await mediabunnyEncoder({ canEncode: p.canEncode }).encode(pcm());
    expect(p.asked.map(([codec]) => codec)).toEqual(['aac', 'opus']);
    expect(isWav(audio)).toBe(true);
  });

  it('the order is an option: Opus first', async () => {
    const p = probe({});
    await mediabunnyEncoder({ canEncode: p.canEncode, prefer: ['opus', 'aac'] }).encode(pcm());
    expect(p.asked.map(([codec]) => codec)).toEqual(['opus', 'aac']);
    const only = probe({});
    await mediabunnyEncoder({ canEncode: only.canEncode, prefer: ['opus'] }).encode(pcm());
    expect(only.asked.map(([codec]) => codec)).toEqual(['opus']);
  });
});
