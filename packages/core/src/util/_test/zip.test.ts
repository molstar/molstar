/**
 * Copyright (c) 2020-2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gzipSync } from 'node:zlib';
import { describe, expect, it, vi } from 'vitest';
import { deflate, inflate, ungzip, unzip, zip, type UnzipLimits } from '../zip/zip.js';
import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';

describe('zip', () => {
  it('roundtrip deflate/inflate', async () => {
    const data = new Uint8Array([1, 2, 3, 4, 5, 6, 7]);
    const deflated = await deflate(SyncRuntimeContext, data);
    const inflated = await inflate(SyncRuntimeContext, deflated);
    expect(inflated).toEqual(data);
  });

  it('roundtrip zip/unzip', async () => {
    const data = {
      'test.foo': new Uint8Array([1, 2, 3, 4, 5, 6, 7]),
    };
    const zipped = await zip(SyncRuntimeContext, data);
    const unzipped = await unzip(SyncRuntimeContext, zipped);
    expect(unzipped).toEqual(data);
  });

  it('ungzip skips an FEXTRA field (a BGZF block) by its XLEN and its bytes', async () => {
    const data = new TextEncoder().encode('data_1ABC\n_entry.id 1ABC\n'.repeat(20));
    const plain = gzipSync(data);
    // BGZF's extra subfield: 'BC', length 2, the block size minus 1
    const extra = [66, 67, 2, 0, 0, 0];
    const withExtra = new Uint8Array(plain.length + 2 + extra.length);
    withExtra.set(plain.subarray(0, 10));
    withExtra[3] |= 4; // FLG.FEXTRA
    withExtra.set([extra.length, 0, ...extra], 10);
    withExtra.set(plain.subarray(10), 10 + 2 + extra.length);
    expect(await ungzip(SyncRuntimeContext, new Uint8Array(plain))).toEqual(data);
    expect(await ungzip(SyncRuntimeContext, withExtra)).toEqual(data);
  });

  // the platform's gzip decoder checks the CRC-32 and the size, with no JS CRC pass over the output (perf-audit-2 01/P1)
  it('ungzip decodes with the platform gzip decoder', async () => {
    const formats: string[] = [];
    const Native = globalThis.DecompressionStream;
    vi.stubGlobal(
      'DecompressionStream',
      class extends Native {
        constructor(format: CompressionFormat) {
          formats.push(format);
          super(format);
        }
      },
    );
    const data = new TextEncoder().encode('data_1ABC\n_entry.id 1ABC\n'.repeat(1000));
    const gz = new Uint8Array(gzipSync(data));
    try {
      expect(await ungzip(SyncRuntimeContext, gz)).toEqual(data);
      expect(formats).toEqual(['gzip']);
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it('ungzip reads a file the platform decoder refuses as before: a wrong CRC-32 is logged', async () => {
    const data = new TextEncoder().encode('data_1ABC\n_entry.id 1ABC\n'.repeat(100));
    const gz = new Uint8Array(gzipSync(data));
    gz[gz.length - 8] ^= 0xff; // the CRC-32's low byte
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    expect(await ungzip(SyncRuntimeContext, gz)).toEqual(data);
    expect(error).toHaveBeenCalledWith("ungzip: checksums don't match");
  });

  describe('untrusted archives', () => {
    const enc = (t: string) => new TextEncoder().encode(t) as Uint8Array<ArrayBuffer>;
    const make = async (extra: Record<string, Uint8Array<ArrayBuffer>> = {}) =>
      new Uint8Array(await zip(SyncRuntimeContext, { 'a.txt': enc('hello hello hello'), ...extra }));
    const read = (z: Uint8Array, limits?: UnzipLimits) =>
      unzip(SyncRuntimeContext, z.buffer as ArrayBuffer, false, limits);

    /** The offsets of an archive's structures: the end record, each central entry and each local header. */
    function layout(z: Uint8Array) {
      const v = new DataView(z.buffer, z.byteOffset, z.byteLength);
      let eocd = z.length - 22;
      while (v.getUint32(eocd, true) !== 0x06054b50) eocd--;
      const central: number[] = [];
      for (let k = 0, o = v.getUint32(eocd + 16, true); k < v.getUint16(eocd + 10, true); k++) {
        central.push(o);
        o += 46 + v.getUint16(o + 28, true) + v.getUint16(o + 30, true) + v.getUint16(o + 32, true);
      }
      return { v, eocd, central, local: central.map((c) => v.getUint32(c + 42, true)) };
    }

    it('reads a valid archive', async () => {
      const r = await read(await make());
      expect(new TextDecoder().decode(r['a.txt'] as Uint8Array)).toBe('hello hello hello');
    });

    it('refuses a truncated archive with no end record', async () => {
      const z = await make();
      await expect(read(z.slice(0, z.length - 22))).rejects.toThrow(/no end of central directory/);
    });

    it('refuses past the limits before inflating: entries and declared sizes', async () => {
      const z = await make({ 'big.bin': new Uint8Array(4 << 20) }); // a few KB deflated
      expect(z.length).toBeLessThan(64 << 10);
      await expect(read(z, { maxEntries: 10, maxBytes: 1 << 20 })).rejects.toThrow(/over the bound of 1048576 bytes/);
      await expect(read(z, { maxEntries: 1, maxBytes: 1 << 30 })).rejects.toThrow(/2 entries, over the bound of 1/);
      expect(Object.keys(await read(z, { maxEntries: 2, maxBytes: 1 << 30 }))).toEqual(['a.txt', 'big.bin']);
    });

    it('refuses a local entry listed twice', async () => {
      const z = await make({ 'b.txt': enc('b') });
      const { v, central, local } = layout(z);
      v.setUint32(central[1] + 42, local[0], true);
      await expect(read(z)).rejects.toThrow(/lists a local entry another lists too/);
    });

    it('refuses an entry that is not its declared size or fails its CRC-32', async () => {
      const stored = await make({ 'a.png': enc('0123456789') }); // .png is stored, not deflated
      const longer = stored.slice();
      const l = layout(longer);
      l.v.setUint32(l.central[1] + 24, 20, true);
      await expect(read(longer)).rejects.toThrow(/runs past|holds 10 bytes, declared 20|fails/);

      const flipped = stored.slice();
      const f = layout(flipped);
      flipped[f.local[1] + 30 + f.v.getUint16(f.local[1] + 26, true)] ^= 1;
      await expect(read(flipped)).rejects.toThrow(/"a\.png" fails its CRC-32/);

      const deflated = await make({ 'big.bin': new Uint8Array(1 << 16) });
      const d = layout(deflated);
      d.v.setUint32(d.central[1] + 24, 1000, true);
      await expect(read(deflated)).rejects.toThrow(/inflates past its declared 1000 bytes|holds/);
    });

    it('reads sizes from the central directory: a local header with zeros (a data descriptor)', async () => {
      const z = await make();
      const { v, local } = layout(z);
      v.setUint16(local[0] + 6, v.getUint16(local[0] + 6, true) | 8, true);
      for (const at of [14, 18, 22]) v.setUint32(local[0] + at, 0, true);
      const r = await read(z);
      expect(new TextDecoder().decode(r['a.txt'] as Uint8Array)).toBe('hello hello hello');
    });

    it('returns what it inflated when the stream is shorter than the buffer', async () => {
      const data = enc('hello hello hello');
      const deflated = await deflate(SyncRuntimeContext, data);
      const inflated = await inflate(SyncRuntimeContext, deflated, new Uint8Array(100));
      expect(inflated.length).toBe(data.length);
    });

    it('reads only the entries `only` accepts, not inflating the rest', async () => {
      const z = await make({ 'state.json': enc('{}') });
      // A damaged entry it skips does not fail the read (its bytes are never looked at)
      const { v, central } = layout(z);
      v.setUint32(central[0] + 16, 0, true); // a.txt's CRC-32
      await expect(read(z)).rejects.toThrow(/fails its CRC-32/);
      const r = await unzip(SyncRuntimeContext, z.buffer as ArrayBuffer, false, undefined, (n) => n === 'state.json');
      expect(Object.keys(r)).toEqual(['state.json']);
      expect(new TextDecoder().decode(r['state.json'] as Uint8Array)).toBe('{}');
    });
  });
});
