/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { describe, expect, it } from 'vitest';
import {
  geminiNarration,
  mp3Duration,
  NarrationError,
  narrationProviders,
  openaiNarration,
  type ProviderAudio,
  pcmToWav,
  refusal,
  wavSamples,
} from '../index.js';
import { base64, mp3Frames, tone, wav } from './helpers.js';

function pcmOf(audio: ProviderAudio): { pcm: Int16Array; sampleRate: number } {
  if (audio.kind !== 'pcm') throw new Error('expected PCM');
  return audio;
}

describe('MP3 length', () => {
  it('is read from the frames (no MP3 is encoded here)', () => {
    expect(Math.abs(mp3Duration(mp3Frames(1.5)) - 1500)).toBeLessThan(30);
    expect(mp3Duration(new Uint8Array(0))).toBe(0);
  });
});

describe('WAV', () => {
  it('only its data chunk is sound: its header and a text chunk at the end are not (they hissed)', () => {
    const pcm = tone(24000, 0.5);
    const w = wavSamples(wav(pcm, 24000));
    expect(w?.sampleRate).toBe(24000);
    expect(w?.pcm.length).toBe(pcm.length);
    expect(w?.pcm[0]).toBe(pcm[0]);
    expect(w?.pcm.at(-1)).toBe(pcm.at(-1));
    // Stereo: mixed to one.
    const stereo = new Int16Array(pcm.length * 2);
    for (let i = 0; i < pcm.length; i++) stereo[i * 2] = stereo[i * 2 + 1] = pcm[i] ?? 0;
    expect(wavSamples(wav(stereo, 24000, 2))?.pcm.length).toBe(pcm.length);
    expect(wavSamples(new Uint8Array([1, 2, 3, 4]))).toBeUndefined();
  });

  it('a WAV written reads back as it was', () => {
    const pcm = tone(16000, 0.25);
    const back = wavSamples(pcmToWav(pcm, 16000));
    expect(back?.sampleRate).toBe(16000);
    expect(Array.from(back?.pcm ?? [])).toEqual(Array.from(pcm));
  });
});

describe('providers', () => {
  it('Gemini: a WAV back (as 3.8 gives it), its sound only as PCM', async () => {
    const pcm = tone(24000, 1);
    const provider = geminiNarration({
      apiKey: 'k',
      fetch: async () =>
        Response.json({
          candidates: [
            { content: { parts: [{ inlineData: { mimeType: 'audio/wav', data: base64(wav(pcm, 24000)) } }] } },
          ],
        }),
    });
    const audio = pcmOf(await provider.synthesize({ text: 'x', voice: 'Kore' }));
    expect(audio.sampleRate).toBe(24000);
    expect(audio.pcm.length).toBe(pcm.length);
    expect(audio.pcm[100]).toBe(pcm[100]);
  });

  it('Gemini: the voice and the delivery note as it takes them; raw PCM back at its rate', async () => {
    const sent: { url: string; body: unknown; key: string | null }[] = [];
    const pcm = tone(16000, 1);
    const bytes = new Uint8Array(pcm.buffer);
    const provider = geminiNarration({
      apiKey: 'test-key',
      fetch: async (url, init) => {
        sent.push({
          url,
          body: JSON.parse(String(init?.body)),
          key: new Headers(init?.headers).get('x-goog-api-key'),
        });
        return Response.json({
          candidates: [
            {
              content: { parts: [{ inlineData: { mimeType: 'audio/L16;codec=pcm;rate=16000', data: base64(bytes) } }] },
            },
          ],
        });
      },
    });
    const audio = pcmOf(
      await provider.synthesize({ text: 'Myoglobin stores oxygen.', voice: 'Charon', style: 'calm' }),
    );
    expect(sent[0]?.url).toContain('models/gemini-3.8-flash-tts:generateContent');
    expect(sent[0]?.key).toBe('test-key');
    expect(sent[0]?.body).toMatchObject({
      contents: [{ parts: [{ text: 'calm:\nMyoglobin stores oxygen.' }] }],
      generationConfig: {
        responseModalities: ['AUDIO'],
        speechConfig: { voiceConfig: { prebuiltVoiceConfig: { voiceName: 'Charon' } } },
      },
    });
    expect(audio.sampleRate).toBe(16000);
    expect(Array.from(audio.pcm)).toEqual(Array.from(pcm));
  });

  it('Gemini: raw PCM with no rate is 24 kHz; a WAV it cannot read, or no audio, is an error in words', async () => {
    const reply = (part: unknown) => async () => Response.json({ candidates: [{ content: { parts: [part] } }] });
    const raw = geminiNarration({
      apiKey: 'k',
      fetch: reply({ inlineData: { mimeType: 'audio/L16', data: base64(new Uint8Array(480)) } }),
    });
    expect(pcmOf(await raw.synthesize({ text: 'x', voice: 'Kore' })).sampleRate).toBe(24000);
    const odd = geminiNarration({
      apiKey: 'k',
      fetch: reply({ inlineData: { mimeType: 'audio/wav', data: base64(new Uint8Array(480)) } }),
    });
    await expect(odd.synthesize({ text: 'x', voice: 'Kore' })).rejects.toThrow(/cannot read/);
    const none = geminiNarration({ apiKey: 'k', fetch: reply({ text: 'no' }) });
    await expect(none.synthesize({ text: 'x', voice: 'Kore' })).rejects.toThrow(/gave no audio/);
  });

  it('OpenAI: MP3 back as it is (encoded, audio/mpeg), its length from its frames; the delivery note as instructions', async () => {
    const mp3 = mp3Frames(2);
    let body: Record<string, unknown> = {};
    const provider = openaiNarration({
      apiKey: 'test-key',
      fetch: async (_url, init) => {
        body = JSON.parse(String(init?.body));
        return new Response(mp3 as Uint8Array<ArrayBuffer>, { headers: { 'content-type': 'audio/mpeg' } });
      },
    });
    const audio = await provider.synthesize({ text: 'Hello.', voice: 'sage', style: 'warm' });
    expect(body).toMatchObject({
      model: 'gpt-4o-mini-tts',
      voice: 'sage',
      input: 'Hello.',
      instructions: 'warm',
      response_format: 'mp3',
    });
    if (audio.kind !== 'encoded') throw new Error('expected encoded audio');
    expect(audio.mimeType).toBe('audio/mpeg');
    expect(audio.bytes).toEqual(mp3);
    expect(Math.abs(audio.durationMs - 2000)).toBeLessThan(40);
  });

  it('a refusal in words, with what to do', async () => {
    const provider = geminiNarration({
      apiKey: 'bad',
      fetch: async () => Response.json({ error: { message: 'API key not valid.' } }, { status: 400 }),
    });
    await expect(provider.synthesize({ text: 'x', voice: 'Kore' })).rejects.toThrow(NarrationError);
    const denied = openaiNarration({
      apiKey: 'bad',
      fetch: async () => Response.json({ error: { message: 'Incorrect API key' } }, { status: 401 }),
    });
    await expect(denied.synthesize({ text: 'x', voice: 'sage' })).rejects.toThrow(
      /401: Incorrect API key\. Check the key\./,
    );
  });

  it('a rate limit passes, an account out of quota does not; no body is just the status', async () => {
    const said = async (status: number, body?: unknown) =>
      (await refusal('X', body === undefined ? new Response('', { status }) : Response.json(body, { status }))).message;
    expect(await said(429, { error: { message: 'Too many requests' } })).toBe(
      'X answered 429: Too many requests. Try again in a moment.',
    );
    expect(await said(429, { error: 'You exceeded your current quota' })).toBe(
      'X answered 429: You exceeded your current quota.',
    );
    expect(await said(500)).toBe('X answered 500.');
  });

  it('the providers there are keys for, in the order offered', () => {
    expect(narrationProviders({}).map((p) => p.id)).toEqual([]);
    expect(narrationProviders({ openai: 'a' }).map((p) => p.id)).toEqual(['openai']);
    expect(narrationProviders({ openai: 'a', google: 'b' }).map((p) => p.id)).toEqual(['google', 'openai']);
  });

  it('the fetch given is the one used', async () => {
    let called = '';
    const [provider] = narrationProviders({ openai: 'k' }, async (url) => {
      called = url;
      return new Response(mp3Frames(0.5) as Uint8Array<ArrayBuffer>);
    });
    await provider?.synthesize({ text: 'x', voice: 'sage' });
    expect(called).toBe('https://api.openai.com/v1/audio/speech');
  });
});
