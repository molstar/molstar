/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

/**
 * Quad items (sphere impostors, points: six vertex ids per item, two triangles, no elements) drawn indexed: four
 * vertex shader invocations per item instead of six (`sphere-impostors.md`, indexed quads). Item `i` is the indices
 * `6i + QuadIndexPattern`, ids of the six-vertex quad whose corners the shaders map (0 (-1, 1), 1 (-1, -1), 2 (1, 1),
 * 4 (1, -1)): its triangles are those of the six vertices, so the shaders, the draw counts (6 per item, indices now)
 * and the draw bases stay, and an indexed draw of a range draws what the unindexed draw of it does
 */
export const QuadIndexPattern = [0, 1, 2, 1, 4, 2] as const;

/** Indices per item, as vertex ids per item unindexed */
export const QuadIndicesPerItem = 6;

const BytesPerItem = QuadIndicesPerItem * 4;
/** Items of the first buffer, so that small scenes do not grow it item by item */
const MinItemCount = 1 << 12;

/**
 * A new array of the indices of the first `itemCount` items: built for an upload and dropped after it, not kept (the
 * buffer holds them; a kept array was as large as the buffer, 59 MB of heap for 2.4 M spheres)
 */
export function createQuadIndexArray(itemCount: number): Uint32Array {
  const indices = new Uint32Array(itemCount * QuadIndicesPerItem);
  for (let o = 0, il = indices.length; o < il; o += QuadIndicesPerItem) {
    for (let k = 0; k < QuadIndicesPerItem; ++k) indices[o + k] = o + QuadIndexPattern[k];
  }
  return indices;
}

export interface QuadIndexBuffer<B> {
  readonly buffer: B;
  /** Changes when the buffer is replaced (grown): draws recorded with an older one are not current */
  readonly id: number;
  readonly itemCount: number;
}

/**
 * The shared index buffer of a context's quad items: one for all, grown by doubling to the largest range drawn, never
 * per item (the pattern is the same for every item, merged or not)
 */
export interface QuadIndices<B> {
  /** A buffer of at least `itemCount` items (`first + count` of a draw), undefined beyond the buffer limit */
  get(itemCount: number): QuadIndexBuffer<B> | undefined;
  /** The id of the current buffer (0 for none): recorded draws binding another are not current */
  readonly id: number;
  readonly byteCount: number;
  dispose(): void;
}

let nextId = 1;

let quadIndexing = true;
/**
 * Whether render items created from now on draw their quads indexed (default); off for comparisons (benchmarks): they
 * draw them unindexed, as past the buffer limit
 */
export function setQuadIndexing(enabled: boolean) {
  quadIndexing = enabled;
}
export function isQuadIndexing() {
  return quadIndexing;
}

/** The most items a buffer of `maxBytes` holds (`GraphicsLimits.maxQuadItemCount`) */
export function getMaxQuadItemCount(maxBytes: number) {
  return Math.floor(maxBytes / BytesPerItem);
}

/**
 * `create` makes a buffer of the indices given, `release` drops one (deferred while recorded work may use it);
 * `maxBytes` the largest buffer the backend makes
 */
export function createQuadIndices<B>(
  create: (indices: Uint32Array) => B,
  release: (buffer: B) => void,
  maxBytes: number,
): QuadIndices<B> {
  const maxItemCount = getMaxQuadItemCount(maxBytes);
  let current: QuadIndexBuffer<B> | undefined;

  return {
    get(itemCount) {
      if (itemCount > maxItemCount) return undefined;
      if (current && current.itemCount >= itemCount) return current;
      const grown = Math.max(itemCount, MinItemCount, current ? current.itemCount * 2 : 0);
      const next = Math.min(grown, maxItemCount);
      if (current) release(current.buffer);
      current = { buffer: create(createQuadIndexArray(next)), id: nextId++, itemCount: next };
      return current;
    },
    get id() {
      return current ? current.id : 0;
    },
    get byteCount() {
      return current ? current.itemCount * BytesPerItem : 0;
    },
    dispose() {
      if (current) release(current.buffer);
      current = undefined;
    },
  };
}
