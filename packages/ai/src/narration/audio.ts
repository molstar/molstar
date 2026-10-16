/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

/**
 * Reading audio as providers send it, and writing the one container with no dependency (WAV). Nothing here encodes
 * MP3: a provider's MP3 is kept as it arrives, and only its length is read (.v6/ai/plan.md).
 */

/**
 * The samples of a WAV file (RIFF): its format (rate, channels, 16-bit) from the "fmt " chunk and
 * only its "data" chunk as sound; other chunks (a LIST of text at the end) are not. Several
 * channels are mixed to one. Undefined when it is not a 16-bit PCM WAV.
 */
export function wavSamples(bytes: Uint8Array): { pcm: Int16Array; sampleRate: number } | undefined {
  const text = (at: number) => String.fromCharCode(...bytes.subarray(at, at + 4));
  if (bytes.length < 12 || text(0) !== 'RIFF' || text(8) !== 'WAVE') return undefined;
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let channels = 1;
  let rate = 0;
  let bits = 0;
  let data: Uint8Array | undefined;
  for (let at = 12; at + 8 <= bytes.length; ) {
    const id = text(at);
    const size = view.getUint32(at + 4, true);
    const body = at + 8;
    if (id === 'fmt ') {
      channels = view.getUint16(body + 2, true);
      rate = view.getUint32(body + 4, true);
      bits = view.getUint16(body + 14, true);
    } else if (id === 'data') {
      // A streamed WAV may give 0 or more than there is: to the end, then.
      const end = size === 0 || body + size > bytes.length ? bytes.length : body + size;
      data = bytes.subarray(body, end);
    }
    at = body + size + (size % 2);
    if (size === 0 && id === 'data') break;
  }
  if (!data || bits !== 16 || !rate || channels < 1) return undefined;
  const all = pcmFromBytes(data);
  if (channels === 1) return { pcm: all, sampleRate: rate };
  const mono = new Int16Array(Math.floor(all.length / channels));
  for (let i = 0; i < mono.length; i++) {
    let sum = 0;
    for (let c = 0; c < channels; c++) sum += all[i * channels + c] ?? 0;
    mono[i] = Math.round(sum / channels);
  }
  return { pcm: mono, sampleRate: rate };
}

/** Little-endian 16-bit samples from raw bytes (as Gemini's audio/L16 arrives). */
export function pcmFromBytes(bytes: Uint8Array): Int16Array {
  const n = Math.floor(bytes.length / 2);
  const view = new DataView(bytes.buffer, bytes.byteOffset, n * 2);
  const out = new Int16Array(n);
  for (let i = 0; i < n; i++) out[i] = view.getInt16(i * 2, true);
  return out;
}

/** Mono 16-bit samples as a WAV file: a 44-byte header (RIFF, "fmt ", "data") and the samples, little-endian. */
export function pcmToWav(pcm: Int16Array, sampleRate: number): Uint8Array {
  const size = pcm.length * 2;
  const out = new Uint8Array(44 + size);
  const view = new DataView(out.buffer);
  const tag = (at: number, text: string) => {
    for (let i = 0; i < 4; i++) out[at + i] = text.charCodeAt(i);
  };
  tag(0, 'RIFF');
  view.setUint32(4, 36 + size, true);
  tag(8, 'WAVE');
  tag(12, 'fmt ');
  view.setUint32(16, 16, true);
  view.setUint16(20, 1, true); // PCM
  view.setUint16(22, 1, true); // mono
  view.setUint32(24, sampleRate, true);
  view.setUint32(28, sampleRate * 2, true);
  view.setUint16(32, 2, true);
  view.setUint16(34, 16, true);
  tag(36, 'data');
  view.setUint32(40, size, true);
  // Written sample by sample: the platform's byte order is not assumed.
  for (let i = 0; i < pcm.length; i++) view.setInt16(44 + i * 2, pcm[i] ?? 0, true);
  return out;
}

/** The length of mono samples at a rate, in ms. */
export function pcmDuration(pcm: Int16Array, sampleRate: number): number {
  return Math.round((pcm.length / sampleRate) * 1000);
}

const BITRATES: Record<string, number[]> = {
  // MPEG-1 Layer III, MPEG-2/2.5 Layer III (kbps by index; 0 and 15 are not valid).
  '1': [0, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320],
  '2': [0, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160],
};
const RATES: Record<number, number[]> = { 3: [44100, 48000, 32000], 2: [22050, 24000, 16000], 0: [11025, 12000, 8000] };

/** An MP3's length, from its frames (Layer III; an ID3 tag skipped), in ms. */
export function mp3Duration(bytes: Uint8Array): number {
  let i = 0;
  if (bytes[0] === 0x49 && bytes[1] === 0x44 && bytes[2] === 0x33) {
    // ID3v2: a 10-byte header, then a syncsafe size.
    const size =
      (((bytes[6] ?? 0) & 0x7f) << 21) |
      (((bytes[7] ?? 0) & 0x7f) << 14) |
      (((bytes[8] ?? 0) & 0x7f) << 7) |
      ((bytes[9] ?? 0) & 0x7f);
    i = 10 + size;
  }
  let seconds = 0;
  while (i + 4 <= bytes.length) {
    const b1 = bytes[i + 1] ?? 0;
    const b2 = bytes[i + 2] ?? 0;
    if (bytes[i] !== 0xff || (b1 & 0xe0) !== 0xe0) {
      i++;
      continue;
    }
    const version = (b1 >> 3) & 0x3; // 3: MPEG-1, 2: MPEG-2, 0: MPEG-2.5
    const layer = (b1 >> 1) & 0x3; // 1: Layer III
    const bitrate = (BITRATES[version === 3 ? '1' : '2'] ?? [])[(b2 >> 4) & 0xf] ?? 0;
    const rate = (RATES[version] ?? [])[(b2 >> 2) & 0x3] ?? 0;
    if (layer !== 1 || version === 1 || !bitrate || !rate) {
      i++;
      continue;
    }
    const padding = (b2 >> 1) & 0x1;
    const samples = version === 3 ? 1152 : 576;
    const length = Math.floor(((samples / 8) * bitrate * 1000) / rate) + padding;
    seconds += samples / rate;
    i += Math.max(1, length);
  }
  return Math.round(seconds * 1000);
}
