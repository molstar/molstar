/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { AssetRangeReader } from '@molstar/core/util/assets-range';
import { type BgzfLimits, BgzfRangeReader, type VirtualOffset } from './internal/bgzf-blocks.js';

// BGZF (`_spec/bgzf.md`): blocked gzip read whole or by virtual offset over core's range reads, each block inflated
// by core's inflaters. The implementation is `internal/bgzf-blocks.ts`.

export {
  BgzfError,
  type BgzfLimits,
  DefaultBgzfLimits,
  inflateBgzf,
  isBgzf,
  VirtualOffset,
} from './internal/bgzf-blocks.js';

export interface BgzfReader {
  /** Inflated bytes from `begin` to `end` (the end block read for its first `within(end)` bytes). */
  read(begin: VirtualOffset, end: VirtualOffset, signal?: AbortSignal): Promise<Uint8Array>;
  /** Whether the file ends with the 28-byte EOF block (one read of its last bytes; false: perhaps cut short). */
  hasEof(signal?: AbortSignal): Promise<boolean>;
  /** Lets go of the inflated blocks it keeps; `ranges` stays the caller's to dispose. */
  dispose(): void;
}

/** A BGZF file read by virtual offset over `ranges` (core's pages), keeping `limits.cacheBlocks` inflated blocks. */
export function bgzfReader(ranges: AssetRangeReader, limits?: Partial<BgzfLimits>): BgzfReader {
  return new BgzfRangeReader(ranges, limits);
}
