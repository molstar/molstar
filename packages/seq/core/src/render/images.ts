/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { LRUCache } from '@molstar/core/util/lru-cache';
import type { ImageData } from './draw-lists.js';

// The shared image cache (_spec/README.md "Render", The shared reduction): every `image` tile a layer draws (a grid's,
// the rows', a navigator's) in one LRU on core's `LRUCache`, keyed by (layer, level, tile, theme version), its capacity
// `imageBytes` over a tile's bytes (tiles are one size). Why one budget: two caches could not lend each other room.

/** The side of an image tile, in samples: 256 × 256. */
export const ImageTile = 256;

/** An image tile's bytes: RGBA a sample. */
export const ImageTileBytes = ImageTile * ImageTile * 4;

interface Slot {
  width: number;
  height: number;
  data: Uint8Array;
  /** The bytes it was made with (`data` is their first `width × height × 4`). */
  bytes: Uint8Array;
}

/**
 * The cache: `get` a tile by its key, else `slot` a tile's image to build into (an evicted tile's bytes when one was let
 * go: no new buffer once the cache is full) and `put` it under its key. A tile is drawn by reference: an image `put`
 * stays valid until it is evicted, which `put` alone does (never while a frame draws it: the frame that drew it laid
 * out before the next `put`).
 */
export class ImageCache {
  private readonly _lru: LRUCache<Slot>;
  private _spare: Slot | undefined;
  private _made = 0;

  /** At most `bytes` over a tile's bytes, at least one tile. */
  constructor(bytes: number, tileBytes = ImageTileBytes) {
    this._lru = LRUCache.create<Slot>(Math.max(1, Math.floor(bytes / tileBytes)));
  }

  /** How many tiles it holds at most. */
  get capacity(): number {
    return this._lru.capacity;
  }

  /** How many tiles it holds. */
  get size(): number {
    return this._lru.entries.count;
  }

  /** How many images it made (the rest reused an evicted tile's bytes): for tests and the allocation tier. */
  get made(): number {
    return this._made;
  }

  get(key: string): ImageData | undefined {
    return LRUCache.get(this._lru, key);
  }

  /** An image of `width` × `height` to build into: the last evicted tile's bytes when they fit, else new ones. */
  slot(width: number, height: number): ImageData {
    const need = width * height * 4;
    const s = this._spare;
    this._spare = undefined;
    if (s && s.bytes.length >= need) {
      s.width = width;
      s.height = height;
      if (s.data.length !== need) s.data = s.bytes.subarray(0, need);
      return s;
    }
    this._made++;
    const bytes = new Uint8Array(Math.max(need, ImageTileBytes));
    const slot: Slot = { width, height, data: bytes.subarray(0, need), bytes };
    return slot;
  }

  /** Keeps `image` under `key` (made by `slot`); the tile it evicts lends its bytes to the next `slot`. */
  put(key: string, image: ImageData): void {
    const out = LRUCache.set(this._lru, key, image as Slot);
    if (out && out !== image) this._spare = out;
  }

  /** Lets every tile go (a theme's tokens changed past their keys, a layer gone). */
  clear(): void {
    this._lru.entries = LRUCache.create<Slot>(this._lru.capacity).entries;
    this._spare = undefined;
  }
}
