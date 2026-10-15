/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { mp3Duration } from './audio.js';
import { type Fetch, type NarrationProvider, type NarrationVoice, refusal } from './provider.js';

const VOICES: NarrationVoice[] = [
  'alloy',
  'ash',
  'ballad',
  'cedar',
  'coral',
  'echo',
  'fable',
  'marin',
  'nova',
  'onyx',
  'sage',
  'shimmer',
  'verse',
].map((id) => ({ id }));

/** OpenAI's text-to-speech (audio/speech): MP3 back, kept as it is; the delivery note as instructions. */
export function openaiNarration(options: { apiKey: string; fetch?: Fetch }): NarrationProvider {
  const doFetch: Fetch = options.fetch ?? ((u, i) => fetch(u, i));
  return {
    id: 'openai',
    name: 'OpenAI',
    // The model ids move to a model catalog later.
    models: [{ id: 'gpt-4o-mini-tts', name: 'GPT-4o mini TTS' }],
    voices: VOICES,
    defaultVoice: 'sage',
    async synthesize({ text, voice, model, style, signal }) {
      const res = await doFetch('https://api.openai.com/v1/audio/speech', {
        method: 'POST',
        headers: { 'content-type': 'application/json', authorization: `Bearer ${options.apiKey}` },
        body: JSON.stringify({
          model: model ?? 'gpt-4o-mini-tts',
          voice,
          input: text,
          response_format: 'mp3',
          ...(style?.trim() ? { instructions: style.trim() } : {}),
        }),
        ...(signal ? { signal } : {}),
      });
      if (!res.ok) throw await refusal('OpenAI', res);
      const bytes = new Uint8Array(await res.arrayBuffer());
      return { kind: 'encoded', bytes, mimeType: 'audio/mpeg', durationMs: mp3Duration(bytes) };
    },
  };
}
