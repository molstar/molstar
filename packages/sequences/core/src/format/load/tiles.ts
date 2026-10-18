/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Tile addressing by region and level (seqspec:sources#the-viewers-side, seqspec:sources#levels-and-the-zoom-gate):
// what the stream contract fixes for an alignment's regions and a 1-D layer's bins, whatever protocol a provider
// speaks. The `seqstar-tiles/1` client is `@molstar/seq-core/stream/tiles`'s (ST2).

/** Samples per tile side; a 1-D tile has as many. */
export const TILE = 256;

/** Addressing along one axis: a level covers 2^level positions per sample, `TILE` samples per tile. No allocation. */
export const Tiles = {
  /** The tile holding `index` at `level`. */
  tileOf: (index: number, level: number) => Math.floor(index / (TILE * 2 ** level)),
  /** The first index a tile covers. */
  start: (tile: number, level: number) => tile * TILE * 2 ** level,
  /** The sample within its tile holding `index`. */
  sampleOf: (index: number, level: number) => Math.floor(index / 2 ** level) % TILE,
  /** Tiles along an axis of `n` positions. */
  count: (n: number, level: number) => Math.ceil(n / (TILE * 2 ** level)),
  /** The finest level at which `n` positions fit in `samples` samples: one tile with `samples` = TILE. */
  levelFor(n: number, samples = TILE) {
    let level = 0;
    while (Math.ceil(n / 2 ** level) > samples) level++;
    return level;
  },
};
