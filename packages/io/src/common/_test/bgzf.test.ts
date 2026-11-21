/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { PCG } from '@molstar/core/data/util/hash-functions';
import { itWithin } from '@molstar/test-harness/timing';
import { BgzfError, bgzfReader, inflateBgzf, isBgzf, VirtualOffset } from '../bgzf.js';
import { blockHeader, walkBlocks } from '../internal/bgzf-blocks.js';
import { bgzip, memoryRanges } from '@molstar/test-harness/genome/bgzf-writer';

const text = (n: number) => {
  const rng = new PCG(7);
  let s = '';
  for (let i = 0; i < n; i++) s += `line ${i}\t${Math.floor(rng.float() * 1e6)}\n`;
  return new TextEncoder().encode(s);
};

/** Only a `BgzfError`, or nothing, comes out of `f`. */
async function noThrowPastAnIssue(f: () => Promise<unknown>) {
  try {
    await f();
  } catch (e) {
    expect((e as Error).name, `${e}`).toBe('BgzfError');
  }
}

describe('VirtualOffset', () => {
  it('packs a block offset and an offset inside it, past 2^32', () => {
    const v = VirtualOffset.of(5_000_000_000, 65535);
    expect(VirtualOffset.block(v)).toBe(5_000_000_000);
    expect(VirtualOffset.within(v)).toBe(65535);
    const view = new DataView(new ArrayBuffer(8));
    view.setBigUint64(0, (5_000_000_000n << 16n) | 65535n, true);
    expect(VirtualOffset.read(view, 0)).toBe(v);
  });

  it('refuses a uint64 past 2^53', () => {
    const view = new DataView(new ArrayBuffer(8));
    view.setBigUint64(0, 1n << 60n, true);
    expect(() => VirtualOffset.read(view, 0)).toThrow(BgzfError);
  });
});

describe('BGZF read whole', () => {
  const data = text(5000);
  const file = bgzip(data, 4000);

  it('tells BGZF from plain gzip by its BC field', () => {
    expect(isBgzf(file.bytes)).toBe(true);
    expect(isBgzf(new Uint8Array(gzipSync(data)))).toBe(false);
    expect(isBgzf(data)).toBe(false);
    expect(isBgzf(file.bytes.subarray(0, 17))).toBe(false);
  });

  it('inflates every block, joined; walks the blocks by BSIZE', async () => {
    expect(await inflateBgzf(file.bytes)).toEqual(data);
    const { offsets, headers } = walkBlocks(file.bytes);
    expect(offsets.slice(0, -1)).toEqual(file.blockOffsets);
    expect(headers.at(-1)!.isize).toBe(0); // the EOF block
    expect(headers.slice(0, -1).reduce((n, h) => n + h.isize, 0)).toBe(data.length);
  });

  it('refuses past maxInflatedBytes before allocating (the ISIZEs summed)', async () => {
    await expect(inflateBgzf(file.bytes, { maxInflatedBytes: data.length - 1 })).rejects.toMatchObject({
      name: 'BgzfError',
      kind: 'too-large',
    });
  });

  it('refuses a file cut short, plain gzip, a wrong CRC-32 and a block declaring past 64 KiB as malformed', async () => {
    const malformed = { name: 'BgzfError', kind: 'malformed' };
    await expect(inflateBgzf(file.bytes.subarray(0, file.bytes.length - 40))).rejects.toMatchObject(malformed);
    await expect(inflateBgzf(new Uint8Array(gzipSync(data)))).rejects.toMatchObject(malformed);
    const flipped = file.bytes.slice();
    flipped[file.blockOffsets[1] - 8] ^= 1; // block 0's CRC-32
    await expect(inflateBgzf(flipped)).rejects.toMatchObject(malformed);
    const big = file.bytes.slice();
    new DataView(big.buffer).setUint32(file.blockOffsets[1] - 4, 70000, true);
    expect(() => blockHeader(big, 0)).toThrow(/declares 70000/);
  });

  it('fuzzed and truncated files never throw past a BgzfError', async () => {
    const rng = new PCG(11);
    const small = bgzip(text(300), 1000).bytes;
    for (let k = 0; k < 200; k++) {
      const f = small.slice(0, k % 3 === 0 ? Math.floor(rng.float() * small.length) : small.length);
      for (let j = 0; j < 1 + (k % 4); j++) f[Math.floor(rng.float() * f.length)] = Math.floor(rng.float() * 256);
      await noThrowPastAnIssue(() => inflateBgzf(f));
      await noThrowPastAnIssue(() => bgzfReader(memoryRanges(f)).read(0, VirtualOffset.of(f.length - 1, 10)));
    }
  });
});

describe('BGZF read by virtual offset', () => {
  const data = text(20000);
  const file = bgzip(data, 3000);

  it('reads across block edges: every range equal to the text between its offsets', async () => {
    const reader = bgzfReader(memoryRanges(file.bytes), { cacheBlocks: 4 });
    const rng = new PCG(3);
    for (let k = 0; k < 100; k++) {
      const a = Math.floor(rng.float() * data.length);
      const b = Math.min(data.length, a + Math.floor(rng.float() * 20000));
      const got = await reader.read(file.voffset(a), file.voffset(b));
      expect(got.length).toBe(b - a);
      expect(got).toEqual(data.subarray(a, b));
    }
    // an end at its block's last byte (within = ISIZE) reads as the next block's start does
    const end = VirtualOffset.of(file.blockOffsets[2], 3000);
    expect(await reader.read(file.voffset(100), end)).toEqual(data.subarray(100, 9000));
    expect(await reader.read(file.voffset(0), file.voffset(data.length))).toEqual(data);
  });

  it('gives the caller its own bytes, and reads nothing for a range whose blocks are cached', async () => {
    const ranges = memoryRanges(file.bytes);
    const reader = bgzfReader(ranges, { cacheBlocks: 8 });
    const a = file.voffset(100);
    const b = file.voffset(900); // inside block 0
    const first = await reader.read(a, b);
    first.fill(0);
    expect(await reader.read(a, b)).toEqual(data.subarray(100, 900));
    const reads = ranges.reads;
    expect(await reader.read(file.voffset(10), file.voffset(7000))).toEqual(data.subarray(10, 7000));
    expect(ranges.reads).toBe(reads + 1); // blocks 1 and 2 were not cached
    expect(await reader.read(file.voffset(20), file.voffset(6000))).toEqual(data.subarray(20, 6000));
    expect(ranges.reads).toBe(reads + 1);
  });

  it('tells the EOF block, and its absence', async () => {
    expect(await bgzfReader(memoryRanges(file.bytes)).hasEof()).toBe(true);
    expect(await bgzfReader(memoryRanges(bgzip(data, 3000, false).bytes)).hasEof()).toBe(false);
  });

  it('refuses a virtual offset off the blocks as malformed', async () => {
    const reader = bgzfReader(memoryRanges(file.bytes));
    await expect(reader.read(0, VirtualOffset.of(file.blockOffsets[1] + 3, 1))).rejects.toMatchObject({
      kind: 'malformed',
    });
    await expect(reader.read(VirtualOffset.of(file.blockOffsets[1], 5000), file.voffset(9000))).rejects.toMatchObject({
      kind: 'malformed',
    });
  });

  it('refuses a read past maxInflatedBytes before a byte is read', async () => {
    const ranges = memoryRanges(file.bytes);
    const reader = bgzfReader(ranges, { maxInflatedBytes: 10000 });
    await expect(reader.read(0, VirtualOffset.of(file.blockOffsets.at(-1)!, 0))).rejects.toMatchObject({
      kind: 'too-large',
    });
    expect(ranges.reads).toBe(0);
    // within the bound
    expect(await reader.read(file.voffset(0), file.voffset(9000))).toEqual(data.subarray(0, 9000));
  });

  it('stops at an aborted signal', async () => {
    const reader = bgzfReader(memoryRanges(file.bytes));
    const ac = new AbortController();
    ac.abort();
    await expect(reader.read(file.voffset(0), file.voffset(50000), ac.signal)).rejects.toThrow();
  });

  // genome-formats.md, Consequences (Performance): a BGZF region of 64 blocks inflated, quick
  const wide = text(80000);
  const wideFile = bgzip(wide, 16384);
  itWithin('a region of 64 blocks inflated', async (within) => {
    const a = 1000;
    const b = a + 64 * 16384 - 2000;
    const got = await within(100, () =>
      bgzfReader(memoryRanges(wideFile.bytes)).read(wideFile.voffset(a), wideFile.voffset(b)),
    );
    expect(got.length).toBe(b - a);
  });
});
