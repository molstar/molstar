/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * The audio encoded behind a port (narration.md): a provider's audio in, the one file the app stores out. The WAV
 * encoder here has no dependency and runs anywhere (Node too); the browser's own AAC and Opus encoders are at
 * `./mediabunny`, so this entry does not pull Mediabunny in.
 */

import { pcmDuration, pcmToWav } from './audio.js';
import type { ProviderAudio } from './provider.js';

/** A narration as it is stored. */
export interface NarrationAudio {
  bytes: Uint8Array;
  /** "audio/wav", "audio/mp4", "audio/webm", "audio/mpeg". */
  mimeType: string;
  /** Without the dot ("wav", "m4a", "webm", "mp3"). */
  extension: string;
  durationMs: number;
}

export interface NarrationEncoder {
  /** Audio already encoded passes through unchanged. */
  encode(audio: ProviderAudio, signal?: AbortSignal): Promise<NarrationAudio>;
}

const EXTENSIONS: Record<string, string> = {
  'audio/mpeg': 'mp3',
  'audio/mp3': 'mp3',
  'audio/mp4': 'm4a',
  'audio/x-m4a': 'm4a',
  'audio/aac': 'aac',
  'audio/webm': 'webm',
  'audio/ogg': 'ogg',
  'audio/opus': 'opus',
  'audio/flac': 'flac',
  'audio/wav': 'wav',
  'audio/wave': 'wav',
  'audio/x-wav': 'wav',
};

/** The file extension for a MIME type (parameters ignored; the subtype when it is not a known one). */
export function extensionOf(mimeType: string): string {
  const type = (mimeType.split(';')[0] ?? '').trim().toLowerCase();
  const known = EXTENSIONS[type];
  if (known) return known;
  const subtype = (type.split('/')[1] ?? '').replace(/^x-/, '').replace(/[^a-z0-9]/g, '');
  return subtype || 'bin';
}

/** A provider's already-encoded audio as it is. */
export function passThrough(audio: Extract<ProviderAudio, { kind: 'encoded' }>): NarrationAudio {
  return {
    bytes: audio.bytes,
    mimeType: audio.mimeType,
    extension: extensionOf(audio.mimeType),
    durationMs: audio.durationMs,
  };
}

/** Mono 16-bit WAV: the larger file (a quarter of what AAC makes), the one that needs nothing. */
export function wavEncoder(): NarrationEncoder {
  return {
    async encode(audio, signal) {
      signal?.throwIfAborted();
      if (audio.kind === 'encoded') return passThrough(audio);
      return {
        bytes: pcmToWav(audio.pcm, audio.sampleRate),
        mimeType: 'audio/wav',
        extension: 'wav',
        durationMs: pcmDuration(audio.pcm, audio.sampleRate),
      };
    },
  };
}
