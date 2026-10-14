/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// PNG for embedded images (seqspec:viewer#the-render-backend): signature, IHDR, IDAT, IEND, CRC-32 and Adler-32, RGBA
// 8-bit, filter 0. Without a `deflate` port the zlib stream is stored (uncompressed) blocks, written straight into the
// one buffer the file is (the pixels are copied once, the checksums computed on the way); with it (the host passes
// `CompressionStream("deflate")` or `zlib.deflateSync`, so this entry imports neither) the scanlines are compressed.
// A size over `MAX_PNG_BYTES` of scanlines is refused before anything is allocated. Platform-free.

/** The most scanline bytes (h · (1 + 4w)) the encoder takes: 512 MiB, an image of about 11,000 × 11,000 pixels. */
export const MAX_PNG_BYTES = 2 ** 29;

const table = new Uint32Array(256);
for (let n = 0; n < 256; n++) {
  let c = n;
  for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
  table[n] = c;
}
/** CRC-32 of `b[from, to)` (the PNG chunk checksum, over the type and data). */
export function crc32(b: Uint8Array, from: number, to: number): number {
  let c = 0xffffffff;
  for (let i = from; i < to; i++) c = table[(c ^ b[i]) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}
const u32 = (b: Uint8Array, at: number, v: number) => {
  b[at] = v >>> 24;
  b[at + 1] = v >>> 16;
  b[at + 2] = v >>> 8;
  b[at + 3] = v;
};
const SIGNATURE = [0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a];
const IEND = [0, 0, 0, 0, 0x49, 0x45, 0x4e, 0x44, 0xae, 0x42, 0x60, 0x82];
const HEAD = 8 + 25 + 8;

/** The signature, IHDR (with its CRC) and IDAT's length and type, in `png`; the IDAT data starts at `HEAD`. */
function header(png: Uint8Array, w: number, h: number, zlen: number): void {
  png.set(SIGNATURE);
  u32(png, 8, 13);
  png.set([0x49, 0x48, 0x44, 0x52], 12);
  u32(png, 16, w);
  u32(png, 20, h);
  png.set([8, 6, 0, 0, 0], 24);
  u32(png, 29, crc32(png, 12, 29));
  u32(png, 33, zlen);
  png.set([0x49, 0x44, 0x41, 0x54], 37);
}
/** The chunk's checksum (over `IDAT` and its data) and the IEND chunk, once the data is in. */
function trailer(png: Uint8Array, zlen: number): Uint8Array {
  u32(png, HEAD + zlen, crc32(png, 37, HEAD + zlen));
  png.set(IEND, HEAD + zlen + 4);
  return png;
}
function check(w: number, h: number, rgba: Uint32Array): number {
  if (!(Number.isInteger(w) && Number.isInteger(h) && w > 0 && h > 0 && w <= 2 ** 31 - 1 && h <= 2 ** 31 - 1))
    throw new Error(`encodePng: ${w} × ${h} is not an image size.`);
  const raw = h * (1 + 4 * w);
  if (raw > MAX_PNG_BYTES) throw new Error(`encodePng: ${w} × ${h} is over the ${MAX_PNG_BYTES / 2 ** 20} MiB bound.`);
  if (rgba.length < w * h) throw new Error(`encodePng: ${rgba.length} pixels for ${w} × ${h}.`);
  return raw;
}

/** The five bytes that open a stored block of `n` bytes. */
function open(png: Uint8Array, at: number, n: number, last: boolean): void {
  png[at] = last ? 1 : 0;
  png[at + 1] = n & 255;
  png[at + 2] = n >>> 8;
  png[at + 3] = ~n & 255;
  png[at + 4] = (~n >>> 8) & 255;
}
/** The PNG with the scanlines stored: one allocation, the pixels copied once, Adler-32 and CRC-32 on the way. */
function stored(w: number, h: number, rgba: Uint32Array, raw: number): Uint8Array {
  const px = new Uint8Array(rgba.buffer, rgba.byteOffset, w * h * 4);
  const blocks = Math.max(1, Math.ceil(raw / 65535));
  const zlen = 2 + raw + 5 * blocks + 4;
  const png = new Uint8Array(HEAD + zlen + 4 + 12);
  header(png, w, h, zlen);
  let at = HEAD;
  png[at++] = 0x78;
  png[at++] = 0x01;
  let a = 1;
  let s = 0;
  let left = 0; // bytes the open block can still take
  let written = 0; // scanline bytes in
  for (let y = 0; y < h; y++) {
    if (left === 0) {
      left = Math.min(65535, raw - written);
      open(png, at, left, written + left === raw);
      at += 5;
    }
    at++; // the filter byte, 0 (the buffer is zeroed)
    left--;
    written++;
    s += a;
    for (let from = y * 4 * w, n = 4 * w; n > 0; ) {
      if (left === 0) {
        left = Math.min(65535, raw - written);
        open(png, at, left, written + left === raw);
        at += 5;
      }
      const k = Math.min(left, n);
      png.set(px.subarray(from, from + k), at);
      for (let i = at, end = at + k; i < end; i++) {
        a += png[i];
        s += a;
      }
      // A block is at most 65,535 bytes, so the sums stay under 2^53 between moduli.
      a %= 65521;
      s %= 65521;
      at += k;
      from += k;
      left -= k;
      written += k;
      n -= k;
    }
  }
  u32(png, at, ((s << 16) | a) >>> 0);
  return trailer(png, zlen);
}

/** The scanlines, each led by filter byte 0: `rgba` is in ImageData's byte order (R, G, B, A in memory). */
function scanlines(w: number, h: number, rgba: Uint32Array, raw: number): Uint8Array {
  const px = new Uint8Array(rgba.buffer, rgba.byteOffset, w * h * 4);
  const out = new Uint8Array(raw);
  for (let y = 0; y < h; y++) out.set(px.subarray(y * 4 * w, (y + 1) * 4 * w), y * (1 + 4 * w) + 1);
  return out;
}
/** The PNG around a zlib stream the host made. */
function wrapped(w: number, h: number, z: Uint8Array): Uint8Array {
  // CMF 8 (deflate) and a header whose 16 bits are a multiple of 31; a raw deflate stream (`"deflate-raw"`) is not.
  if (z.length < 6 || (z[0] & 15) !== 8 || ((z[0] << 8) | z[1]) % 31 !== 0)
    throw new Error('encodePng: `deflate` must return a zlib stream (CompressionStream("deflate"), zlib.deflateSync).');
  const png = new Uint8Array(HEAD + z.length + 4 + 12);
  header(png, w, h, z.length);
  png.set(z, HEAD);
  return trailer(png, z.length);
}

/** A PNG of w × h pixels; with `deflate` (a zlib stream out) compressed, else stored. */
export function encodePng(
  w: number,
  h: number,
  rgba: Uint32Array,
  deflate?: (bytes: Uint8Array) => Uint8Array | Promise<Uint8Array>,
): Uint8Array | Promise<Uint8Array> {
  const raw = check(w, h, rgba);
  if (!deflate) return stored(w, h, rgba, raw);
  const z = deflate(scanlines(w, h, rgba, raw));
  return z instanceof Uint8Array ? wrapped(w, h, z) : z.then((s) => wrapped(w, h, s));
}

const B64 = Uint8Array.from("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/", (c) => c.charCodeAt(0));
/** Standard base64 of `b` (the platform's `btoa` and `Buffer` are not in this entry). */
export function base64(b: Uint8Array): string {
  const out = new Uint8Array(Math.ceil(b.length / 3) * 4).fill(61);
  let o = 0;
  for (let i = 0; i < b.length; i += 3) {
    const n = (b[i] << 16) | ((b[i + 1] ?? 0) << 8) | (b[i + 2] ?? 0);
    out[o++] = B64[n >>> 18];
    out[o++] = B64[(n >>> 12) & 63];
    if (i + 1 < b.length) out[o++] = B64[(n >>> 6) & 63];
    else o++;
    if (i + 2 < b.length) out[o++] = B64[n & 63];
    else o++;
  }
  return new TextDecoder().decode(out);
}
