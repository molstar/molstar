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

/**
 * Box items (cylinder impostors: six vertex ids per item, the corners of the box the shader maps from `id % 6`), four
 * triangles: item `i` is the indices `6i + BoxIndexPattern`, 12 per item (the draw counts and bases count indices)
 */
export const BoxIndexPattern = [0, 1, 2, 1, 4, 2, 2, 4, 3, 4, 5, 3] as const;

/**
 * Line items (wide lines: four vertex ids per item, the corners left and right of the start and the end; text: the
 * corners of a glyph, background or tether): two triangles, item `i` the indices `4i + LineIndexPattern`
 */
export const LineIndexPattern = [0, 1, 2, 1, 3, 2] as const;

/** Items drawn by a shared index pattern: `pattern` the indices of an item, of its `verticesPerItem` vertex ids */
export interface ItemIndexPattern {
  readonly name: 'quad' | 'box' | 'line';
  readonly pattern: readonly number[];
  readonly verticesPerItem: number;
}

export const QuadItems: ItemIndexPattern = { name: 'quad', pattern: QuadIndexPattern, verticesPerItem: 6 };
export const BoxItems: ItemIndexPattern = { name: 'box', pattern: BoxIndexPattern, verticesPerItem: 6 };
export const LineItems: ItemIndexPattern = { name: 'line', pattern: LineIndexPattern, verticesPerItem: 4 };

/** Items of the first buffer, so that small scenes do not grow it item by item */
const MinItemCount = 1 << 12;

/**
 * A new array of the indices of the first `itemCount` items of `items`: built for an upload and dropped after it, not
 * kept (the buffer holds them; a kept array was as large as the buffer, 59 MB of heap for 2.4 M spheres)
 */
export function createItemIndexArray(itemCount: number, items: ItemIndexPattern): Uint32Array {
  const { pattern, verticesPerItem } = items;
  const perItem = pattern.length;
  const indices = new Uint32Array(itemCount * perItem);
  for (let i = 0, o = 0; i < itemCount; ++i, o += perItem) {
    const v = i * verticesPerItem;
    for (let k = 0; k < perItem; ++k) indices[o + k] = v + pattern[k];
  }
  return indices;
}

/** The indices of the first `itemCount` quads (`createItemIndexArray`) */
export function createQuadIndexArray(itemCount: number): Uint32Array {
  return createItemIndexArray(itemCount, QuadItems);
}

export interface QuadIndexBuffer<B> {
  readonly buffer: B;
  /** Changes when the buffer is replaced (grown): draws recorded with an older one are not current */
  readonly id: number;
  readonly itemCount: number;
}

/**
 * The shared index buffer of a context's items of one pattern (quads, boxes): one for all, grown by doubling to the
 * largest range drawn, never per item (the pattern is the same for every item, merged or not)
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

/**
 * Whether items of `items` draw indexed by the shared indices of their pattern: boxes and lines always (they have no
 * unindexed form), quads unless turned off (`setQuadIndexing`)
 */
export function isItemIndexing(items: ItemIndexPattern | undefined): items is ItemIndexPattern {
  return !!items && (items.name !== 'quad' || quadIndexing);
}

/** The shared indices of `items` of a context's resources */
export function getItemIndices<B>(
  resources: {
    readonly quadIndices: QuadIndices<B>;
    readonly boxIndices: QuadIndices<B>;
    readonly lineIndices: QuadIndices<B>;
  },
  items: ItemIndexPattern,
): QuadIndices<B> {
  if (items.name === 'box') return resources.boxIndices;
  if (items.name === 'line') return resources.lineIndices;
  return resources.quadIndices;
}

/** The most items of `items` a draw of the shared indices reaches (`GraphicsLimits`) */
export function getMaxItemCount(
  limits: { readonly maxQuadItemCount: number; readonly maxBoxItemCount: number; readonly maxLineItemCount: number },
  items: ItemIndexPattern,
) {
  if (items.name === 'box') return limits.maxBoxItemCount;
  if (items.name === 'line') return limits.maxLineItemCount;
  return limits.maxQuadItemCount;
}

/** The most items of `items` a buffer of `maxBytes` holds (`GraphicsLimits.maxQuadItemCount`, `maxBoxItemCount`) */
export function getMaxQuadItemCount(maxBytes: number, items: ItemIndexPattern = QuadItems) {
  return Math.floor(maxBytes / (items.pattern.length * 4));
}

/**
 * `create` makes a buffer of the indices given, `release` drops one (deferred while recorded work may use it);
 * `maxBytes` the largest buffer the backend makes; `items` the pattern (quads by default)
 */
export function createQuadIndices<B>(
  create: (indices: Uint32Array) => B,
  release: (buffer: B) => void,
  maxBytes: number,
  items: ItemIndexPattern = QuadItems,
): QuadIndices<B> {
  const maxItemCount = getMaxQuadItemCount(maxBytes, items);
  const bytesPerItem = items.pattern.length * 4;
  let current: QuadIndexBuffer<B> | undefined;

  return {
    get(itemCount) {
      if (itemCount > maxItemCount) return undefined;
      if (current && current.itemCount >= itemCount) return current;
      const grown = Math.max(itemCount, MinItemCount, current ? current.itemCount * 2 : 0);
      const next = Math.min(grown, maxItemCount);
      if (current) release(current.buffer);
      current = { buffer: create(createItemIndexArray(next, items)), id: nextId++, itemCount: next };
      return current;
    },
    get id() {
      return current ? current.id : 0;
    },
    get byteCount() {
      return current ? current.itemCount * bytesPerItem : 0;
    },
    dispose() {
      if (current) release(current.buffer);
      current = undefined;
    },
  };
}
