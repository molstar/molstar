/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/** A 440 Hz tone, 16-bit at `rate`. */
export function tone(rate: number, seconds = 1): Int16Array {
  const out = new Int16Array(Math.round(rate * seconds));
  for (let i = 0; i < out.length; i++) out[i] = Math.round(8000 * Math.sin((2 * Math.PI * 440 * i) / rate));
  return out;
}

export const base64 = (bytes: Uint8Array) => {
  let s = '';
  for (const b of bytes) s += String.fromCharCode(b);
  return btoa(s);
};

/** A WAV as Gemini 3.8 gives it: RIFF, "fmt ", "data", then a LIST chunk of text at the end. */
export function wav(pcm: Int16Array, rate: number, channels = 1): Uint8Array {
  const data = new Uint8Array(pcm.buffer, pcm.byteOffset, pcm.byteLength);
  const list = new TextEncoder().encode('INFOISFT\u000e\u0000\u0000\u0000Lavf61.7.100\u0000\u0000');
  const out = new Uint8Array(12 + 24 + 8 + data.length + 8 + list.length);
  const v = new DataView(out.buffer);
  const tag = (at: number, t: string) => {
    for (let i = 0; i < 4; i++) out[at + i] = t.charCodeAt(i);
  };
  tag(0, 'RIFF');
  v.setUint32(4, out.length - 8, true);
  tag(8, 'WAVE');
  tag(12, 'fmt ');
  v.setUint32(16, 16, true);
  v.setUint16(20, 1, true);
  v.setUint16(22, channels, true);
  v.setUint32(24, rate, true);
  v.setUint32(28, rate * 2 * channels, true);
  v.setUint16(32, 2 * channels, true);
  v.setUint16(34, 16, true);
  tag(36, 'data');
  v.setUint32(40, data.length, true);
  out.set(data, 44);
  tag(44 + data.length, 'LIST');
  v.setUint32(48 + data.length, list.length, true);
  out.set(list, 52 + data.length);
  return out;
}

/**
 * MPEG-1 Layer III frames, silent: `seconds` of them at 44.1 kHz, 128 kbps (417 bytes a frame, 1152 samples), which is
 * all `mp3Duration` reads. A stand-in for a provider's MP3 (nothing here encodes one).
 */
export function mp3Frames(seconds: number): Uint8Array {
  const frames = Math.round((seconds * 44100) / 1152);
  const out = new Uint8Array(frames * 417);
  for (let f = 0; f < frames; f++) out.set([0xff, 0xfb, 0x90, 0x00], f * 417);
  return out;
}

export const text = (bytes: Uint8Array, at: number, n = 4) => String.fromCharCode(...bytes.subarray(at, at + n));
