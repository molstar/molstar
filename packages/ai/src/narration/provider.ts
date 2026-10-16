/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Voicing a scene's narration (narration.md): a text-to-speech provider behind one interface, so none is forced
 * (Gemini is one; OpenAI another; more can join). A provider returns the audio as it comes: raw PCM, or already
 * encoded (a provider's MP3 is kept as it arrives); the encoder (./encoder.ts) turns that into the one kind of file
 * the app stores. No MP3 is ever encoded here.
 */

export interface NarrationVoice {
  id: string;
  /** A word for how it sounds, as the provider says ("Firm", "Warm"). */
  note?: string;
}

export interface NarrationModel {
  id: string;
  name: string;
}

export interface NarrationRequest {
  text: string;
  voice: string;
  /** Default: the provider's first model. */
  model?: string;
  /** How to say it, in words ("calm, documentary"): passed as the provider takes it. */
  style?: string;
  signal?: AbortSignal;
}

/** What a provider gives back, before it is encoded. */
export type ProviderAudio =
  /** Mono 16-bit samples. */
  | { kind: 'pcm'; pcm: Int16Array; sampleRate: number }
  /** A file as the provider made it (an MP3, say), with its length. */
  | { kind: 'encoded'; bytes: Uint8Array; mimeType: string; durationMs: number };

export interface NarrationProvider {
  /** As the app names providers ("google", "openai"): its key is that provider's key. */
  id: string;
  name: string;
  /** The first is the default. */
  models: readonly NarrationModel[];
  voices: readonly NarrationVoice[];
  defaultVoice: string;
  synthesize(request: NarrationRequest): Promise<ProviderAudio>;
}

/** A provider's failure, in words for the person (what the provider said, and what to do). */
export class NarrationError extends Error {
  override name = 'NarrationError';
}

export type Fetch = (url: string, init?: RequestInit) => Promise<Response>;

/** What a provider said when it refused (its own message when it gives one). */
export async function refusal(provider: string, res: Response): Promise<NarrationError> {
  let detail = '';
  try {
    const body = (await res.json()) as { error?: { message?: string } | string };
    detail = typeof body.error === 'string' ? body.error : (body.error?.message ?? '');
  } catch {
    detail = '';
  }
  // A rate limit passes; an account out of credit or quota does not.
  const spent = /credit|quota|billing/i.test(detail);
  const hint =
    res.status === 401 || res.status === 403
      ? ' Check the key.'
      : res.status === 429 && !spent
        ? ' Try again in a moment.'
        : '';
  const said = detail.trim() ? `: ${detail.trim()}${/[.!?]$/.test(detail.trim()) ? '' : '.'}` : '.';
  return new NarrationError(`${provider} answered ${res.status}${said}${hint}`);
}
