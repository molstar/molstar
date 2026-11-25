/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Task } from '@molstar/core/task';
import type { Asset, AssetManager } from './assets.js';
import { wholeReader } from './internal/range-pages.js';

// Range reads of an asset (`_spec/assets.md` §10): byte ranges of a URL or a file, read by page and kept in a page
// cache by asset and page, shared by the readers of the asset and let go with the last of them. Exported through
// `assets.ts`. The implementation (`DefaultAssetManager.openRange`'s page stores, the whole-asset fallback) is in
// `internal/range-pages.ts`, which the package does not publish.

/** A byte range: `[offset, offset + length)`. */
export interface AssetRange {
  readonly offset: number;
  readonly length: number;
}

export interface AssetRangeReader {
  /** The asset's size in bytes, known after open (`Content-Range` total, `stat`, `Blob.size`). */
  readonly size: number;
  /**
   * One result per range, in order, each a copy-free view of cached pages where one page covers it. Ranges past the
   * end are clipped; an empty result is a zero-length array. Rejects as aborted with the signal's reason.
   */
  read(ranges: readonly AssetRange[], signal?: AbortSignal): Promise<Uint8Array[]>;
  /** Gives back the reference: pages are let go with the last holder. */
  dispose(): void;
}

export interface RangeReadOptions
  extends Pick<AssetManager.ResolveOptions, 'credentials' | 'redirect' | 'signal' | 'maxBytes'> {
  /** The page size, bytes (default 64 KiB; a power of two, 4 KiB to 4 MiB). */
  pageSize?: number;
  /** Pages held at most for this asset (default 256: 16 MiB); least recently used go first. */
  maxPages?: number;
}

/**
 * A range read cannot be trusted: a `Content-Range` that is not the one asked for or has no total, a 416 on a valid
 * range, the asset's size, `ETag` or `Last-Modified` changed between pages, a body longer or shorter than its range, a
 * reader used after its manager was cleared. Matched by `name` across copies of the core, as `SizeLimitExceeded` is.
 */
export class AssetRangeError extends Error {
  constructor(
    why: string,
    /** The URL or the file name. */
    readonly source?: string,
  ) {
    super(source ? `${source}: ${why}` : why);
    this.name = 'AssetRangeError';
  }
}

/**
 * `manager.openRange` when the manager has it; else the whole asset is resolved (`binary`, `shared`, bounded by
 * `maxBytes`) and ranges are served as views of its bytes: correct, not lazy.
 */
export function openRanges(
  manager: AssetManager,
  asset: Asset,
  options: RangeReadOptions = {},
): Task<AssetRangeReader> {
  if (manager.openRange) return manager.openRange(asset, options);
  return wholeReader(manager, asset, options);
}
