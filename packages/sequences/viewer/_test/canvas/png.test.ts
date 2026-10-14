/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The PNG writer in Chromium (SV1b; seqspec:viewer#the-render-backend): the browser's own decoders read what it writes.
// Chunk CRCs are recomputed here bit by bit (not the writer's table), `DecompressionStream("deflate")` inflates the
// zlib stream (it checks the Adler-32), and the image decoder gives the pixels back; stored and compressed through
// `CompressionStream`.

import { describe, expect, it } from "vitest";
import { compressionDeflate } from "../../src/dom/index.js";
import { encodePng } from "../../src/svg/index.js";

/** CRC-32, one bit at a time. */
function slowCrc(b: Uint8Array): number {
  let c = 0xffffffff;
  for (const byte of b) {
    c ^= byte;
    for (let k = 0; k < 8; k++) c = c & 1 ? (c >>> 1) ^ 0xedb88320 : c >>> 1;
  }
  return (c ^ 0xffffffff) >>> 0;
}
async function inflate(z: Uint8Array): Promise<Uint8Array> {
  const stream = new Blob([z as Uint8Array<ArrayBuffer>]).stream().pipeThrough(new DecompressionStream("deflate"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}
/** The chunks checked, the zlib stream inflated: [w, h, scanlines]. */
async function read(png: Uint8Array, checkCrc: boolean): Promise<[number, number, Uint8Array]> {
  const dv = new DataView(png.buffer, png.byteOffset, png.length);
  let [at, w, h] = [8, 0, 0];
  const z: Uint8Array[] = [];
  while (at < png.length) {
    const len = dv.getUint32(at);
    const type = String.fromCharCode(...png.subarray(at + 4, at + 8));
    if (checkCrc) expect(dv.getUint32(at + 8 + len), `${type}'s CRC`).toBe(slowCrc(png.subarray(at + 4, at + 8 + len)));
    if (type === "IHDR") [w, h] = [dv.getUint32(at + 8), dv.getUint32(at + 12)];
    if (type === "IDAT") z.push(png.subarray(at + 8, at + 8 + len));
    at += 12 + len;
  }
  return [w, h, await inflate(new Uint8Array(await new Blob(z as BlobPart[]).arrayBuffer()))];
}
/** What the image decoder makes of the file: RGBA bytes. */
async function decoded(png: Uint8Array, w: number, h: number): Promise<Uint8ClampedArray> {
  const bmp = await createImageBitmap(new Blob([png as Uint8Array<ArrayBuffer>], { type: "image/png" }), {
    premultiplyAlpha: "none",
    colorSpaceConversion: "none",
  });
  const c = new OffscreenCanvas(w, h);
  const g = c.getContext("2d", { willReadFrequently: true }) as OffscreenCanvasRenderingContext2D;
  g.drawImage(bmp, 0, 0);
  return g.getImageData(0, 0, w, h).data;
}
/** Opaque or fully transparent pixels in a pattern: exact through any alpha handling. */
function pattern(w: number, h: number): Uint32Array {
  const px = new Uint32Array(w * h);
  for (let y = 0; y < h; y++)
    for (let x = 0; x < w; x++) {
      const clear = (x + y) % 7 === 0;
      // ImageData's byte order, R first in memory: little-endian ABGR.
      px[y * w + x] = clear
        ? 0
        : ((255 << 24) | (((y * 5) & 255) << 16) | (((x * 3) & 255) << 8) | ((x + y) & 255)) >>> 0;
    }
  return px;
}

describe.each([
  ["stored", undefined],
  ["through CompressionStream", compressionDeflate],
])("a PNG written %s", (_, deflate) => {
  it.each([
    [1, 1],
    [37, 19],
    [16384, 3], // a row longer than a stored block
    [200, 700],
  ])("%i × %i: the CRCs, the zlib stream and the image decoder agree with the pixels", async (w, h) => {
    const px = pattern(w, h);
    const png = await encodePng(w, h, px, deflate);
    const [rw, rh, raw] = await read(png, w * h < 100_000);
    expect([rw, rh]).toEqual([w, h]);
    const row = 1 + 4 * w;
    expect(raw.length).toBe(h * row);
    const bytes = new Uint8Array(px.buffer);
    for (let y = 0; y < h; y++) {
      expect(raw[y * row]).toBe(0);
      expect(raw.subarray(y * row + 1, (y + 1) * row)).toEqual(bytes.subarray(y * 4 * w, (y + 1) * 4 * w));
    }
    const img = await decoded(png, w, h);
    expect(new Uint8Array(img.buffer, img.byteOffset, img.length)).toEqual(
      // Transparent pixels come back as black or as written, depending on the decoder: compare the opaque ones' bytes.
      new Uint8Array(img.length).map((_, i) => ((px[i >> 2] >>> 24) & 255 ? bytes[i] : img[i])),
    );
  });
});

describe("a damaged file", () => {
  it("a wrong Adler-32 or a wrong CRC is refused by the browser's decoders", async () => {
    const png = (await encodePng(8, 8, pattern(8, 8))) as Uint8Array;
    const adler = png.slice();
    adler[adler.length - 12 - 4 - 1] ^= 0xff; // the last byte of the zlib stream's Adler-32
    await expect(read(adler, false)).rejects.toBeDefined();
    const crc = png.slice();
    crc[crc.length - 12 - 4] ^= 0xff; // IDAT's CRC
    await expect(decoded(crc, 8, 8)).rejects.toBeDefined();
  });
});

describe("a large image", () => {
  it("2,048 × 2,048 stored: written in bounded time and read by the browser", async () => {
    const [w, h] = [2048, 2048];
    const px = new Uint32Array(w * h).fill(0xff336699);
    const t = performance.now();
    const png = await encodePng(w, h, px);
    expect(performance.now() - t).toBeLessThan(3000); // about 100 ms
    expect(png.length).toBeLessThan(w * h * 4 * 1.001 + 4000);
    const img = await decoded(png, w, h);
    expect([img[0], img[1], img[2], img[3]]).toEqual([0x99, 0x66, 0x33, 0xff]);
    expect(img.slice(4 * (w * h - 1))).toEqual(new Uint8ClampedArray([0x99, 0x66, 0x33, 0xff]));
    // Compressed, the same image is about 21 KB (deflate's 1,032:1 limit).
    expect((await encodePng(w, h, px, compressionDeflate)).length).toBeLessThan(40_000);
  });
});
