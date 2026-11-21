/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { SyncRuntimeContext } from '@molstar/core/task/execution/synchronous';
import type { AssetRangeReader } from '@molstar/core/util/assets-range';
import { LRUCache } from '@molstar/core/util/lru-cache';
import { crc } from '@molstar/core/util/zip/checksum';
import { inflateStream } from '@molstar/core/util/zip/inflate';
import { checkCompressionStreamSupport } from '@molstar/core/util/zip/util';
import { inflateRaw } from '@molstar/core/util/zip/zip';

// BGZF's implementation (`_spec/bgzf.md`): virtual offsets, the block walk (BSIZE), a block inflated with core's
// inflaters, the range reader. Not published (`./common/internal/*`): `common/bgzf.ts` is the surface, and tests reach
// the block walk here.

/** A block's file offset × 2^16 + an offset into its inflated bytes; a number (files under 2^37 bytes: 128 GiB). */
export type VirtualOffset = number;

const Within = 0x10000;

export const VirtualOffset = {
  of(block: number, within: number): VirtualOffset {
    return block * Within + within;
  },
  block(v: VirtualOffset): number {
    return Math.floor(v / Within);
  },
  within(v: VirtualOffset): number {
    return v % Within;
  },
  /** A little-endian uint64 at `at` (a `.tbi`'s), refused past 2^53. */
  read(view: DataView, at: number): VirtualOffset {
    const lo = view.getUint32(at, true);
    const hi = view.getUint32(at + 4, true);
    if (hi >= 0x200000) throw new BgzfError('malformed', `a virtual offset at byte ${at} is past 2^53`);
    return hi * 0x100000000 + lo;
  },
};

export interface BgzfLimits {
  /** Inflated bytes one read gives, its blocks together (default 64 MiB: the store's answer bound). */
  readonly maxInflatedBytes: number;
  /** Inflated blocks a reader keeps for a pan back (default 16). */
  readonly cacheBlocks: number;
}

export const DefaultBgzfLimits: BgzfLimits = { maxInflatedBytes: 64 * 1024 * 1024, cacheBlocks: 16 };

/** Not BGZF, a block's `BC` field, size, CRC-32 or ISIZE wrong, or a bound passed; by `name` across copies. */
export class BgzfError extends Error {
  constructor(
    readonly kind: 'malformed' | 'too-large',
    message: string,
  ) {
    super(message);
    this.name = 'BgzfError';
  }
}

/** A block's bound, compressed and inflated (BGZF: BSIZE is 16 bits, a block inflates to at most 64 KiB). */
export const MaxBlockBytes = 0x10000;

/** The 28-byte empty block bgzip writes last. */
export const EofBlock = new Uint8Array([
  0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0,
]);

export interface BlockHeader {
  /** The block's whole size, compressed: BSIZE + 1. */
  readonly size: number;
  /** Where its deflate stream starts, from the block's start. */
  readonly data: number;
  /** Its inflated size (ISIZE). */
  readonly isize: number;
  /** Its CRC-32. */
  readonly crc32: number;
}

/**
 * The header of the block at `at` in `bytes`; undefined when `bytes` end before its header or its end (`need` says
 * how many bytes from `at` it needs). Throws `malformed` for a block that is not one (SAM spec §4.1).
 */
export function blockHeader(bytes: Uint8Array, at: number): BlockHeader | { readonly need: number } {
  const n = bytes.length - at;
  if (n < 18) return { need: 18 };
  if (bytes[at] !== 0x1f || bytes[at + 1] !== 0x8b || bytes[at + 2] !== 8 || !(bytes[at + 3] & 4)) {
    throw new BgzfError('malformed', `no BGZF block at byte ${at} (not gzip with an extra field)`);
  }
  const xlen = bytes[at + 10] | (bytes[at + 11] << 8);
  if (n < 12 + xlen) return { need: 12 + xlen };
  let bsize = -1;
  for (let s = at + 12, end = at + 12 + xlen; s + 4 <= end; ) {
    const len = bytes[s + 2] | (bytes[s + 3] << 8);
    if (bytes[s] === 66 && bytes[s + 1] === 67 && len === 2 && s + 6 <= end) bsize = bytes[s + 4] | (bytes[s + 5] << 8);
    s += 4 + len;
  }
  if (bsize < 0) throw new BgzfError('malformed', `the gzip member at byte ${at} has no BC field (not BGZF)`);
  const size = bsize + 1;
  const data = 12 + xlen;
  if (size < data + 8) throw new BgzfError('malformed', `the block at byte ${at} is smaller than its header`);
  if (n < size) return { need: size };
  const v = at + size - 8;
  const crc32 = (bytes[v] | (bytes[v + 1] << 8) | (bytes[v + 2] << 16) | (bytes[v + 3] << 24)) >>> 0;
  const isize = (bytes[v + 4] | (bytes[v + 5] << 8) | (bytes[v + 6] << 16) | (bytes[v + 7] << 24)) >>> 0;
  if (isize > MaxBlockBytes)
    throw new BgzfError('malformed', `the block at byte ${at} declares ${isize} bytes inflated`);
  return { size, data, isize, crc32 };
}

/** Whether `bytes` start with a BGZF block's header: gzip, FEXTRA, a `BC` subfield. */
export function isBgzf(bytes: Uint8Array): boolean {
  if (bytes.length < 18 || bytes[0] !== 0x1f || bytes[1] !== 0x8b || bytes[2] !== 8 || !(bytes[3] & 4)) return false;
  const end = Math.min(12 + (bytes[10] | (bytes[11] << 8)), bytes.length);
  for (let s = 12; s + 4 <= end; s += 4 + (bytes[s + 2] | (bytes[s + 3] << 8))) {
    if (bytes[s] === 66 && bytes[s + 1] === 67 && bytes[s + 2] === 2 && bytes[s + 3] === 0) return true;
  }
  return false;
}

/** The blocks of `bytes` from `at` (a whole number of them, else `malformed`: cut short), as `[offset, header]`. */
export function walkBlocks(bytes: Uint8Array, at = 0): { offsets: number[]; headers: BlockHeader[] } {
  const offsets: number[] = [];
  const headers: BlockHeader[] = [];
  while (at < bytes.length) {
    const h = blockHeader(bytes, at);
    if (!('size' in h)) throw new BgzfError('malformed', `the block at byte ${at} is cut short`);
    offsets.push(at);
    headers.push(h);
    at += h.size;
  }
  return { offsets, headers };
}

/**
 * One block inflated, into `into` when given (`header.isize` long): the platform's gzip decoder (which checks the
 * CRC-32 and ISIZE), else core's JS inflater and a CRC-32 here. Never a second inflater.
 */
export async function inflateBlock(
  block: Uint8Array,
  header: BlockHeader,
  into?: Uint8Array<ArrayBuffer>,
): Promise<Uint8Array<ArrayBuffer>> {
  const out = into ?? new Uint8Array(header.isize);
  if (header.isize === 0) return out;
  const member = block.subarray(0, header.size) as Uint8Array<ArrayBuffer>;
  let got: Uint8Array<ArrayBuffer>;
  try {
    if (checkCompressionStreamSupport('gzip')) {
      got = await inflateStream(SyncRuntimeContext, member, 'gzip', out);
    } else {
      got = await inflateRaw(SyncRuntimeContext, member.subarray(header.data, header.size - 8), out);
      if (got.length === header.isize && crc(got, 0, got.length) >>> 0 !== header.crc32) {
        throw new Error('its CRC-32 is wrong');
      }
    }
  } catch (e) {
    throw new BgzfError('malformed', `a BGZF block does not inflate: ${(e as Error)?.message ?? e}`);
  }
  if (got.length !== header.isize) {
    throw new BgzfError('malformed', `a BGZF block inflates to ${got.length} bytes, declared ${header.isize}`);
  }
  return got;
}

export function limitsOf(limits?: Partial<BgzfLimits>): BgzfLimits {
  return { ...DefaultBgzfLimits, ...limits };
}

/** A BGZF file read whole: every block inflated in turn, joined; the sum of ISIZEs bounded before allocating. */
export async function inflateBgzf(
  bytes: Uint8Array,
  limits?: Partial<BgzfLimits>,
  signal?: AbortSignal,
): Promise<Uint8Array> {
  const { maxInflatedBytes } = limitsOf(limits);
  if (!isBgzf(bytes)) throw new BgzfError('malformed', 'not BGZF (no gzip block with a BC field)');
  const { offsets, headers } = walkBlocks(bytes);
  let total = 0;
  for (const h of headers) total += h.isize;
  if (total > maxInflatedBytes) {
    throw new BgzfError('too-large', `inflates to ${total} bytes, over the bound of ${maxInflatedBytes}`);
  }
  const out = new Uint8Array(total);
  for (let i = 0, o = 0; i < headers.length; o += headers[i].isize, i++) {
    signal?.throwIfAborted();
    await inflateBlock(bytes.subarray(offsets[i]), headers[i], out.subarray(o, o + headers[i].isize));
  }
  return out;
}

/** A block read by the range reader: its inflated bytes and where the next block starts. */
interface Block {
  readonly bytes: Uint8Array;
  readonly next: number;
}

/** The range reader of `bgzfReader`; `readWithin` and `blockAt` are for tabix (a per-call bound, the header walk). */
export class BgzfRangeReader {
  private readonly limits: BgzfLimits;
  private cache: LRUCache<Block>;
  private eof: Promise<boolean> | undefined;
  private disposed = false;

  constructor(
    private readonly ranges: AssetRangeReader,
    limits?: Partial<BgzfLimits>,
  ) {
    this.limits = limitsOf(limits);
    this.cache = LRUCache.create(this.limits.cacheBlocks);
  }

  /** The file's size, compressed. */
  get size() {
    return this.ranges.size;
  }

  read(begin: VirtualOffset, end: VirtualOffset, signal?: AbortSignal): Promise<Uint8Array> {
    return this.readWithin(begin, end, this.limits.maxInflatedBytes, signal);
  }

  /** `read` within `maxBytes` inflated (at most the reader's bound): refused `too-large` before inflating. */
  async readWithin(
    begin: VirtualOffset,
    end: VirtualOffset,
    maxBytes: number,
    signal?: AbortSignal,
  ): Promise<Uint8Array> {
    if (this.disposed) throw new BgzfError('malformed', 'the BGZF reader is disposed');
    const bound = Math.min(maxBytes, this.limits.maxInflatedBytes);
    if (!(begin >= 0 && end >= begin)) throw new BgzfError('malformed', `a read from ${begin} to ${end}`);
    const first = VirtualOffset.block(begin);
    const last = VirtualOffset.block(end);
    const startWithin = VirtualOffset.within(begin);
    const lastWithin = VirtualOffset.within(end);
    // the blocks read: from the first to the last, the last only when some of its bytes are
    const stop = lastWithin > 0 ? last : last - 1;
    if (stop < first) return new Uint8Array(0);
    if (first >= this.size) throw new BgzfError('malformed', `a virtual offset past the file's ${this.size} bytes`);
    // compressed bytes are at most about the inflated ones (a block that does not shrink is stored): the span is
    // bounded before it is read, so an index's offsets never size a read
    if (last - first > bound + MaxBlockBytes) {
      throw new BgzfError('too-large', `a read of ${last - first} compressed bytes, over the bound of ${bound}`);
    }
    const [bytes] = await this.ranges.read(
      [{ offset: first, length: Math.min(last + MaxBlockBytes, this.size) - first }],
      signal,
    );
    // the blocks' headers first (cached blocks as they are), the exact length bounded, then inflated
    const blocks: { at: number; isize: number; header?: BlockHeader; cached?: Block }[] = [];
    let at = first;
    while (at <= stop) {
      const cached = LRUCache.get(this.cache, `${at}`);
      if (cached) {
        blocks.push({ at, isize: cached.bytes.length, cached });
        at = cached.next;
        continue;
      }
      const h = blockHeader(bytes, at - first);
      if (!('size' in h)) throw new BgzfError('malformed', `the block at byte ${at} is cut short`);
      blocks.push({ at, isize: h.isize, header: h });
      at += h.size;
    }
    const tail = blocks[blocks.length - 1];
    if (lastWithin > 0 ? tail.at !== last : at !== last) {
      throw new BgzfError('malformed', `no block starts at byte ${last} (a virtual offset off the blocks)`);
    }
    if (startWithin > blocks[0].isize) {
      throw new BgzfError('malformed', `the virtual offset ${begin} is past its block's ${blocks[0].isize} bytes`);
    }
    if (lastWithin > tail.isize) {
      throw new BgzfError('malformed', `the virtual offset ${end} is past its block's ${tail.isize} bytes`);
    }
    let length = -startWithin;
    for (const b of blocks) length += b.isize;
    if (lastWithin > 0) length -= tail.isize - lastWithin;
    if (length > bound) throw new BgzfError('too-large', `a read of ${length} bytes, over the bound of ${bound}`);
    if (length <= 0) return new Uint8Array(0);

    const inflated: Uint8Array[] = [];
    for (const b of blocks) {
      signal?.throwIfAborted();
      if (b.cached) {
        inflated.push(b.cached.bytes);
        continue;
      }
      const got = await inflateBlock(bytes.subarray(b.at - first), b.header!);
      LRUCache.set(this.cache, `${b.at}`, { bytes: got, next: b.at + b.header!.size });
      inflated.push(got);
    }
    if (inflated.length === 1) return inflated[0].subarray(startWithin, startWithin + length);
    const out = new Uint8Array(length);
    for (let i = 0, o = 0; i < inflated.length; i++) {
      const from = i === 0 ? startWithin : 0;
      const to = i === inflated.length - 1 && lastWithin > 0 ? lastWithin : inflated[i].length;
      out.set(inflated[i].subarray(from, to), o);
      o += to - from;
    }
    return out;
  }

  /** The block at file offset `at` (inflated, cached) and where the next starts; undefined at the file's end. */
  async blockAt(at: number, signal?: AbortSignal): Promise<Block | undefined> {
    if (this.disposed) throw new BgzfError('malformed', 'the BGZF reader is disposed');
    if (at >= this.size) return;
    const cached = LRUCache.get(this.cache, `${at}`);
    if (cached) return cached;
    const [bytes] = await this.ranges.read([{ offset: at, length: Math.min(MaxBlockBytes, this.size - at) }], signal);
    const h = blockHeader(bytes, 0);
    if (!('size' in h)) throw new BgzfError('malformed', `the block at byte ${at} is cut short`);
    const block = { bytes: await inflateBlock(bytes, h), next: at + h.size };
    LRUCache.set(this.cache, `${at}`, block);
    return block;
  }

  hasEof(signal?: AbortSignal): Promise<boolean> {
    if (this.disposed) return Promise.reject(new BgzfError('malformed', 'the BGZF reader is disposed'));
    if (!this.eof) {
      const size = this.size;
      this.eof =
        size < EofBlock.length
          ? Promise.resolve(false)
          : this.ranges.read([{ offset: size - EofBlock.length, length: EofBlock.length }], signal).then(([tail]) => {
              for (let i = 0; i < EofBlock.length; i++) if (tail[i] !== EofBlock[i]) return false;
              return true;
            });
      // an aborted check is not kept
      this.eof.catch(() => {
        this.eof = undefined;
      });
    }
    return this.eof;
  }

  /** Lets go of the cached blocks; the range reader stays its owner's (who disposes it). */
  dispose() {
    this.disposed = true;
    this.cache = LRUCache.create(1);
  }
}
