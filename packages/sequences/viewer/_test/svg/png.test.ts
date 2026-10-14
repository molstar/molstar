/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The PNG writer (SV1b; seqspec:viewer#the-render-backend) in Node: chunk CRCs and the zlib stream's Adler-32 checked
// by zlib itself, the pixels back, stored blocks across their 65,535-byte boundaries, a host's `deflate` (sync and
// async) and a host that returns the wrong stream, the size bound, the stored encoder's time on a large image, and
// `crc32` against zlib's.

import { deflateSync, crc32 as zlibCrc32 } from "node:zlib";
import { describe, expect, it } from "vitest";
import { itWithin } from "../../../core/_test/model/timing.js";
import { encodePng, MAX_PNG_BYTES } from "../../src/svg/index.js";
import { crc32 } from "../../src/svg/png.js";
import { decodePng, pixelsOf } from "./png.js";

/** Pixels that differ everywhere (a stored block's boundary lands in a different channel of a different pixel). */
function noise(n: number): Uint32Array {
  const px = new Uint32Array(n);
  let s = 12345;
  for (let i = 0; i < n; i++) {
    s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
    px[i] = s;
  }
  return px;
}
const bytesOf = (px: Uint32Array) => new Uint8Array(px.buffer, px.byteOffset, px.length * 4);
/** The same bytes (a deep match of a large image's took most of its test). */
const same = (a: Uint8Array, b: Uint8Array) => Buffer.from(a.buffer, a.byteOffset, a.length).equals(b);

describe("the checksums", () => {
  // The stored stream's Adler-32 is computed on the way and checked by inflate in every case below.
  it("crc32 equals zlib's", () => {
    for (const n of [0, 1, 2, 255, 5551, 5552, 5553, 65535, 200_000]) {
      const b = bytesOf(noise(Math.ceil(n / 4))).subarray(0, n);
      expect(crc32(b, 0, n)).toBe(zlibCrc32(b));
    }
  });
});

describe("stored blocks", () => {
  // Rows of 4w + 1 bytes against blocks of 65,535: widths that put a boundary mid-pixel, on a row end, and a block per row.
  it.each([
    [1, 1],
    [1, 3],
    [3, 5],
    [16383, 4], // a row is 65,533 bytes
    [16384, 3], // 65,537: a row longer than a block
    [32767, 5],
    [100, 700],
    [1, 40_000],
    [64, 255], // 257-byte rows, 65,535 bytes: exactly one block, final
    [64, 510], // exactly two
    [54, 400], // 217-byte rows: row 302's filter byte is the first block's last byte
    [4, 3856], // 17-byte rows: a block ends on a row end, and one more row follows
  ])("%i × %i: every chunk's CRC and the zlib stream hold, and the pixels come back", async (w, h) => {
    const px = noise(w * h);
    const png = (await encodePng(w, h, px)) as Uint8Array;
    const d = decodePng(png);
    expect([d.w, d.h, d.chunks]).toEqual([w, h, ["IHDR", "IDAT", "IEND"]]);
    expect(same(pixelsOf(d), bytesOf(px))).toBe(true);
    // One zlib stream of stored blocks: the exact size, a block per 65,535 bytes.
    const raw = h * (1 + 4 * w);
    expect(d.zlib.length).toBe(2 + raw + 5 * Math.ceil(raw / 65535) + 4);
    expect(png.length).toBe(8 + 25 + 12 + d.zlib.length + 12);
  });

  it("an image that is a view into a larger buffer", async () => {
    const big = noise(100);
    const view = big.subarray(10, 10 + 6 * 5);
    expect(pixelsOf(decodePng((await encodePng(6, 5, view)) as Uint8Array))).toEqual(bytesOf(view));
  });
});

describe("a host's deflate", () => {
  it("compresses (sync or async) and the result is the same pixels", async () => {
    const px = new Uint32Array(300 * 200).fill(0xff3366aa);
    const stored = (await encodePng(300, 200, px)) as Uint8Array;
    for (const deflate of [(b: Uint8Array) => deflateSync(b), async (b: Uint8Array) => deflateSync(b)]) {
      const png = await encodePng(300, 200, px, deflate);
      expect(same(pixelsOf(decodePng(png)), bytesOf(px))).toBe(true);
      expect(png.length).toBeLessThan(stored.length / 20);
    }
  });
  it("refuses a stream that is not zlib (a raw deflate, nothing)", () => {
    const px = new Uint32Array(4);
    expect(() => encodePng(2, 2, px, () => new Uint8Array([1, 2, 3, 4, 5, 6, 7, 8]))).toThrow(/zlib stream/);
    expect(() => encodePng(2, 2, px, () => new Uint8Array(0))).toThrow(/zlib stream/);
  });
});

describe("the bounds", () => {
  it("refuses a size that is not an image, too few pixels, and anything over the bound before allocating", () => {
    for (const [w, h] of [
      [0, 1],
      [1, 0],
      [-1, 2],
      [1.5, 2],
      [Number.NaN, 2],
    ])
      expect(() => encodePng(w, h, new Uint32Array(4))).toThrow(/not an image size/);
    expect(() => encodePng(3, 3, new Uint32Array(4))).toThrow(/4 pixels for 3 × 3/);
    // 20,000 × 20,000 is 1.6 GB of scanlines: refused on its size alone (the pixels are never looked at).
    expect(() => encodePng(20_000, 20_000, new Uint32Array(1))).toThrow(/over the 512 MiB bound/);
    expect(MAX_PNG_BYTES).toBe(2 ** 29);
  });

  itWithin("a large image is written stored in a time and a single buffer, not copies of it", async (within) => {
    const [w, h] = [3000, 3000]; // 36 MB of pixels
    const px = new Uint32Array(w * h);
    for (let i = 0; i < px.length; i += 11) px[i] = i;
    // About 150 ms; the bound is for a loaded machine.
    const png = (await within(3000, () => encodePng(w, h, px))) as Uint8Array;
    expect(png.buffer.byteLength).toBe(png.length); // the one buffer is the file
    expect(png.length).toBeLessThan(36_000_000 * 1.001 + 3000 * 5);
    const d = decodePng(png);
    expect([d.w, d.h]).toEqual([w, h]);
    expect(pixelsOf(d).subarray(4 * 11, 4 * 12)).toEqual(bytesOf(px).subarray(4 * 11, 4 * 12));
  });
});
