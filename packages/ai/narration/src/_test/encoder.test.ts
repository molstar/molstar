/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  extensionOf,
  type NarrationAudio,
  type NarrationProvider,
  narrate,
  type ProviderAudio,
  passThrough,
  wavEncoder,
  wavSamples,
} from '../index.js';
import { mp3Frames, text, tone } from './helpers.js';

describe('wavEncoder', () => {
  it('PCM is a valid 16-bit mono WAV with its length', async () => {
    const pcm = tone(24000, 1.5);
    const audio = await wavEncoder().encode({ kind: 'pcm', pcm, sampleRate: 24000 });
    expect(audio.mimeType).toBe('audio/wav');
    expect(audio.extension).toBe('wav');
    expect(audio.durationMs).toBe(1500);
    const b = audio.bytes;
    const v = new DataView(b.buffer, b.byteOffset, b.byteLength);
    expect(b.length).toBe(44 + pcm.length * 2);
    expect(text(b, 0)).toBe('RIFF');
    expect(v.getUint32(4, true)).toBe(b.length - 8);
    expect(text(b, 8)).toBe('WAVE');
    expect(text(b, 12)).toBe('fmt ');
    expect(v.getUint16(20, true)).toBe(1);
    expect(v.getUint16(22, true)).toBe(1);
    expect(v.getUint32(24, true)).toBe(24000);
    expect(v.getUint32(28, true)).toBe(48000);
    expect(v.getUint16(32, true)).toBe(2);
    expect(v.getUint16(34, true)).toBe(16);
    expect(text(b, 36)).toBe('data');
    expect(v.getUint32(40, true)).toBe(pcm.length * 2);
    expect(Array.from(wavSamples(b)?.pcm ?? [])).toEqual(Array.from(pcm));
  });

  it('an empty sound is a header and no samples', async () => {
    const audio = await wavEncoder().encode({ kind: 'pcm', pcm: new Int16Array(0), sampleRate: 24000 });
    expect(audio.bytes.length).toBe(44);
    expect(audio.durationMs).toBe(0);
  });

  it('an MP3 passes through unchanged', async () => {
    const bytes = mp3Frames(1);
    const audio = await wavEncoder().encode({ kind: 'encoded', bytes, mimeType: 'audio/mpeg', durationMs: 1000 });
    expect(audio.bytes).toBe(bytes);
    expect(audio).toMatchObject({ mimeType: 'audio/mpeg', extension: 'mp3', durationMs: 1000 });
  });

  it('an abort stops it', async () => {
    const stop = new AbortController();
    stop.abort();
    await expect(
      wavEncoder().encode({ kind: 'pcm', pcm: tone(8000, 0.1), sampleRate: 8000 }, stop.signal),
    ).rejects.toThrow();
  });
});

describe('extensions', () => {
  it('from the MIME type, its parameters ignored', () => {
    expect(extensionOf('audio/mpeg')).toBe('mp3');
    expect(extensionOf('audio/mp4')).toBe('m4a');
    expect(extensionOf('audio/webm;codecs=opus')).toBe('webm');
    expect(extensionOf('Audio/X-WAV')).toBe('wav');
    expect(extensionOf('audio/x-foo')).toBe('foo');
    expect(extensionOf('')).toBe('bin');
    expect(
      passThrough({ kind: 'encoded', bytes: new Uint8Array(1), mimeType: 'audio/ogg', durationMs: 1 }).extension,
    ).toBe('ogg');
  });
});

describe('narrate', () => {
  const provider = (audio: ProviderAudio, seen: unknown[] = []): NarrationProvider => ({
    id: 'stub',
    name: 'Stub',
    models: [{ id: 'm', name: 'M' }],
    voices: [{ id: 'v' }],
    defaultVoice: 'v',
    async synthesize(request) {
      seen.push(request);
      return audio;
    },
  });

  it('synthesizes, then encodes: PCM becomes a WAV', async () => {
    const seen: unknown[] = [];
    const audio = await narrate(provider({ kind: 'pcm', pcm: tone(24000, 1), sampleRate: 24000 }, seen), wavEncoder(), {
      text: 'Hello.',
      voice: 'v',
    });
    expect(seen).toEqual([{ text: 'Hello.', voice: 'v' }]);
    expect(audio).toMatchObject({ mimeType: 'audio/wav', extension: 'wav', durationMs: 1000 });
    expect(text(audio.bytes, 0)).toBe('RIFF');
  });

  it('an encoded sound is kept as it is', async () => {
    const bytes = mp3Frames(1);
    const audio: NarrationAudio = await narrate(
      provider({ kind: 'encoded', bytes, mimeType: 'audio/mpeg', durationMs: 1000 }),
      wavEncoder(),
      { text: 'x', voice: 'v' },
    );
    expect(audio.bytes).toBe(bytes);
    expect(audio.extension).toBe('mp3');
  });

  it("the request's signal reaches the encoder", async () => {
    const stop = new AbortController();
    const encoded: (AbortSignal | undefined)[] = [];
    await narrate(
      provider({ kind: 'pcm', pcm: tone(8000, 0.1), sampleRate: 8000 }),
      {
        async encode(audio, signal) {
          encoded.push(signal);
          return wavEncoder().encode(audio, signal);
        },
      },
      { text: 'x', voice: 'v', signal: stop.signal },
    );
    expect(encoded).toEqual([stop.signal]);
  });

  it('a provider that fails fails the whole', async () => {
    const failing: NarrationProvider = {
      ...provider({ kind: 'pcm', pcm: new Int16Array(0), sampleRate: 8000 }),
      synthesize: async () => {
        throw new Error('no key');
      },
    };
    await expect(narrate(failing, wavEncoder(), { text: 'x', voice: 'v' })).rejects.toThrow('no key');
  });
});
