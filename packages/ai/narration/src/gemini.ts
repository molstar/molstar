/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { pcmFromBytes, wavSamples } from './audio.js';
import { type Fetch, NarrationError, type NarrationProvider, type NarrationVoice, refusal } from './provider.js';

/** Gemini's prebuilt voices, with the word Google gives each. */
const VOICES: NarrationVoice[] = [
  ['Kore', 'Firm'],
  ['Charon', 'Informative'],
  ['Sulafat', 'Warm'],
  ['Achird', 'Friendly'],
  ['Iapetus', 'Clear'],
  ['Erinome', 'Clear'],
  ['Schedar', 'Even'],
  ['Vindemiatrix', 'Gentle'],
  ['Sadaltager', 'Knowledgeable'],
  ['Rasalgethi', 'Informative'],
  ['Gacrux', 'Mature'],
  ['Algieba', 'Smooth'],
  ['Despina', 'Smooth'],
  ['Achernar', 'Soft'],
  ['Zephyr', 'Bright'],
  ['Autonoe', 'Bright'],
  ['Puck', 'Upbeat'],
  ['Laomedeia', 'Upbeat'],
  ['Aoede', 'Breezy'],
  ['Callirrhoe', 'Easy-going'],
  ['Umbriel', 'Easy-going'],
  ['Zubenelgenubi', 'Casual'],
  ['Orus', 'Firm'],
  ['Alnilam', 'Firm'],
  ['Pulcherrima', 'Forward'],
  ['Sadachbia', 'Lively'],
  ['Leda', 'Youthful'],
  ['Fenrir', 'Excitable'],
  ['Enceladus', 'Breathy'],
  ['Algenib', 'Gravelly'],
].map(([id, note]) => ({ id: id as string, note: note as string }));

/**
 * Gemini's text-to-speech (generateContent with audio out): a WAV or raw 16-bit PCM (audio/L16 at its rate) back,
 * given as PCM for the encoder. The delivery note goes in the text, as Gemini takes it.
 */
export function geminiNarration(options: { apiKey: string; fetch?: Fetch }): NarrationProvider {
  const doFetch: Fetch = options.fetch ?? ((u, i) => fetch(u, i));
  return {
    id: 'google',
    name: 'Google Gemini',
    // The model ids move to a model catalog later.
    models: [
      { id: 'gemini-3.8-flash-tts', name: 'Gemini 3.8 Flash TTS' },
      { id: 'gemini-3.8-flash-lite-tts', name: 'Gemini 3.8 Flash Lite TTS' },
    ],
    voices: VOICES,
    defaultVoice: 'Kore',
    async synthesize({ text, voice, model, style, signal }) {
      const id = model ?? 'gemini-3.8-flash-tts';
      const res = await doFetch(`https://generativelanguage.googleapis.com/v1beta/models/${id}:generateContent`, {
        method: 'POST',
        headers: { 'content-type': 'application/json', 'x-goog-api-key': options.apiKey },
        body: JSON.stringify({
          contents: [{ parts: [{ text: style?.trim() ? `${style.trim()}:\n${text}` : text }] }],
          generationConfig: {
            responseModalities: ['AUDIO'],
            speechConfig: { voiceConfig: { prebuiltVoiceConfig: { voiceName: voice } } },
          },
        }),
        ...(signal ? { signal } : {}),
      });
      if (!res.ok) throw await refusal('Gemini', res);
      const body = (await res.json()) as {
        candidates?: { content?: { parts?: { inlineData?: { mimeType?: string; data?: string } }[] } }[];
      };
      const part = body.candidates?.[0]?.content?.parts?.find((p) => p.inlineData?.data)?.inlineData;
      if (!part?.data) throw new NarrationError('Gemini gave no audio for this text.');
      const bytes = Uint8Array.from(atob(part.data), (c) => c.charCodeAt(0));
      // A WAV (as Gemini 3.8 gives it: its header and a text chunk are not sound), or raw 16-bit PCM
      // (audio/L16, as earlier models gave it).
      const wav = wavSamples(bytes);
      if (wav) return { kind: 'pcm', pcm: wav.pcm, sampleRate: wav.sampleRate };
      if (/wav/i.test(part.mimeType ?? '')) throw new NarrationError('Gemini gave audio in a form this cannot read.');
      const rate = Number(/rate=(\d+)/.exec(part.mimeType ?? '')?.[1] ?? 24000);
      return { kind: 'pcm', pcm: pcmFromBytes(bytes), sampleRate: rate };
    },
  };
}
