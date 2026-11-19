/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { adler, crc } from '@molstar/core/util/zip/checksum';

// A PNG of RGBA pixels, synchronous and without a DOM (_spec/README.md "Render", Backends: the SVG's images): the
// rows filtered `none`, the zlib stream in stored blocks. Why stored: the SVG is written synchronously, core's deflate
// is asynchronous, and an image tile is at most 256 × 256 (256 KB raw).

const Signature = [137, 80, 78, 71, 13, 10, 26, 10];

/** A stored block's most bytes. */
const StoredMax = 65535;

function chunk(out: Uint8Array, at: number, type: string, data: Uint8Array): number {
  const v = new DataView(out.buffer, out.byteOffset);
  v.setUint32(at, data.length);
  for (let i = 0; i < 4; i++) out[at + 4 + i] = type.charCodeAt(i);
  out.set(data, at + 8);
  v.setUint32(at + 8 + data.length, crc(out, at + 4, 4 + data.length) >>> 0);
  return at + 12 + data.length;
}

/** `width` × `height` RGBA pixels (`rgba`, row-major) as a PNG file's bytes. */
export function encodePng(width: number, height: number, rgba: Uint8Array): Uint8Array {
  const row = width * 4 + 1;
  const raw = new Uint8Array(row * height);
  for (let y = 0; y < height; y++) raw.set(rgba.subarray(y * width * 4, (y + 1) * width * 4), y * row + 1);
  const blocks = Math.max(1, Math.ceil(raw.length / StoredMax));
  const z = new Uint8Array(2 + raw.length + 5 * blocks + 4);
  z[0] = 0x78;
  z[1] = 0x01;
  let o = 2;
  for (let b = 0; b < blocks; b++) {
    const start = b * StoredMax;
    const n = Math.min(StoredMax, raw.length - start);
    z[o] = b === blocks - 1 ? 1 : 0;
    z[o + 1] = n & 0xff;
    z[o + 2] = n >>> 8;
    z[o + 3] = ~n & 0xff;
    z[o + 4] = (~n >>> 8) & 0xff;
    z.set(raw.subarray(start, start + n), o + 5);
    o += 5 + n;
  }
  new DataView(z.buffer).setUint32(o, adler(raw, 0, raw.length) >>> 0);
  const ihdr = new Uint8Array(13);
  const h = new DataView(ihdr.buffer);
  h.setUint32(0, width);
  h.setUint32(4, height);
  ihdr[8] = 8; // bit depth
  ihdr[9] = 6; // RGBA
  const out = new Uint8Array(8 + 25 + 12 + z.length + 12);
  out.set(Signature, 0);
  let at = chunk(out, 8, 'IHDR', ihdr);
  at = chunk(out, at, 'IDAT', z);
  chunk(out, at, 'IEND', new Uint8Array(0));
  return out;
}

const Base64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';

/** Bytes as base64 (no `btoa`, no `Buffer`: the same in every runtime). */
export function base64Of(bytes: Uint8Array): string {
  const parts: string[] = [];
  let s = '';
  for (let i = 0; i < bytes.length; i += 3) {
    const a = bytes[i];
    const b = i + 1 < bytes.length ? bytes[i + 1] : 0;
    const c = i + 2 < bytes.length ? bytes[i + 2] : 0;
    s +=
      Base64[a >> 2] +
      Base64[((a & 3) << 4) | (b >> 4)] +
      (i + 1 < bytes.length ? Base64[((b & 15) << 2) | (c >> 6)] : '=') +
      (i + 2 < bytes.length ? Base64[c & 63] : '=');
    if (s.length >= 8192) {
      parts.push(s);
      s = '';
    }
  }
  parts.push(s);
  return parts.join('');
}

/** An image's pixels as a PNG data URL. */
export function pngDataUrl(width: number, height: number, rgba: Uint8Array): string {
  return `data:image/png;base64,${base64Of(encodePng(width, height, rgba))}`;
}
