/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { Interval, OrderedSet } from '@molstar/core/data/int';
import { applyMarkerAction, applyMarkerActionAtPosition, MarkerAction } from '@molstar/core/util/marker-action';
import type { GraphicsContext } from '@molstar/graphics/gl/context';
import { CopyMarkerBlock } from '@molstar/graphics/gl/culled-data';
import { CopyMarkerTexels } from '@molstar/graphics/gl/instance-expansion';
import { StorePageSize } from './page-store.js';

// The markers of a page store's copies (a nested list's, `paged.ts`): its template particles drawn once a copy from one
// store index each, its store markers (by store index) every copy's at once. Per copy one marker for all of it, or a
// block of markers by its template's store indices where marked in part (made as it is, let go once nothing in it is
// marked: a hover's copy at a time, not every copy's particles). What the GPU reads (`sCopyMarkers`, the vertex stage
// by an instance's copy, `CopyMarkerBlock`): two words a copy, then the blocks' bytes. The blocks within a budget (the
// device's, `getCopyBlockBudget`): past it a copy marked in part is marked whole (its word, no block), what is drawn
// coarser, the CPU's selection as it was. Why: a selection differing by copy makes a block each (an entity in half the
// isg cell's 4,096 copies 1.1 GB), past the device's limit the buffer not made

/** A template of the store as its copies' markers need it: its pages (store pages) and its copies' count */
export interface CopyMarkerTemplate {
  readonly pages: ArrayLike<number>;
  readonly copyCount: number;
}

export interface CopyMarkers {
  /** the copies' count: the templates' copies in their order (the cull's copy indices) */
  readonly copyCount: number;
  /** the copies with something marked */
  readonly markedCount: number;
  /** the copies marked whole for want of a block (past the budget) since the layout */
  readonly coarseCount: number;
  /** two words a copy (`CopyMarkerBlock`) */
  readonly words: Uint32Array;
  /** the blocks' bytes, after the words in the GPU's buffer (their offsets in the words from its start) */
  readonly bytes: Uint8Array;
  /** the GPU's buffer's bytes */
  readonly byteLength: number;
  /**
   * the blocks' bytes changed since the last call (`[start, end)` pairs in `bytes`, sorted, merged), cleared;
   * undefined: all of them (laid out again)
   */
  takeChanged(): number[] | undefined;
  /** the templates, every marker cleared */
  layout(templates: readonly CopyMarkerTemplate[]): void;
  /** `action` on copy `c` whole */
  markCopy(c: number, action: MarkerAction): boolean;
  /** `action` on the particles of copy `c` at `indices` (its template's store indices) */
  markIndices(c: number, indices: OrderedSet<number>, action: MarkerAction): boolean;
  /** `action` on every copy whole */
  markAll(action: MarkerAction): boolean;
  /**
   * the markers `values` holds for copy `c`'s template's particles (by store index: the store's, a template's marks
   * every copy's) or-ed into its own (a block made where any): a template's marks made each copy's
   */
  include(c: number, values: Uint8Array): boolean;
  /** whether copy `c` has something marked */
  isMarked(c: number): boolean;
  /** the marker of copy `c`'s particle at store index `index` */
  get(c: number, index: number): number;
  /**
   * where copy `c` is marked: undefined where marked whole (no block: its marker every particle's), null where not
   * known (a block filled whole or past a sixteenth marked), else positions in its block (from its template's first
   * store index; some may be clear again)
   */
  getMarkedPositions(c: number): readonly number[] | null | undefined;
}

/** changed bytes closer than this one range (a write each, what lies between less) */
const ChangedGap = 4096;
/** changed ranges handed out at most: more merged across their smallest gaps */
const MaxChangedRanges = 64;

/** `[start, end)` pairs sorted and merged, at most `MaxChangedRanges` (across the smallest gaps) */
function mergeChanged(pairs: number[]) {
  const order = Array.from({ length: pairs.length / 2 }, (_, k) => k).sort((a, b) => pairs[2 * a] - pairs[2 * b]);
  const merged: number[] = [];
  for (const k of order) {
    const start = pairs[2 * k];
    const end = pairs[2 * k + 1];
    const n = merged.length;
    if (n > 0 && start <= merged[n - 1] + ChangedGap) merged[n - 1] = Math.max(merged[n - 1], end);
    else merged.push(start, end);
  }
  while (merged.length > 2 * MaxChangedRanges) {
    let best = 1;
    for (let k = 3; k < merged.length - 1; k += 2) {
      if (merged[k + 1] - merged[k] < merged[best + 1] - merged[best]) best = k;
    }
    merged.splice(best, 2);
  }
  return merged;
}

/** the blocks' bytes without a device (a core WebGPU device's binding limit, 128 MiB, halved) */
export const DefaultCopyBlockBudget = 64 * 1024 * 1024;

/**
 * The blocks' bytes a device takes: half what its copy markers' binding holds (a storage buffer's limit on WebGPU, the
 * copy texture's rows on WebGL: `CopyMarkerTexels` texels of four bytes a row), room for the buffer's growth by half
 * again and the words
 */
export function getCopyBlockBudget(graphics: GraphicsContext | undefined) {
  if (!graphics) return DefaultCopyBlockBudget;
  const limit = graphics.compute
    ? graphics.compute.maxStorageBufferSize
    : CopyMarkerTexels * 4 * graphics.resources.limits.maxTextureSize;
  return Math.floor(limit / 2 / 4) * 4;
}

/** `getBudget` the blocks' bytes at most, read at each layout (a representation's graphics set later) */
export function createCopyMarkers(getBudget: () => number = () => DefaultCopyBlockBudget): CopyMarkers {
  let budget = getBudget();
  let copyCount = 0;
  let markedCount = 0;
  let coarseCount = 0;
  /** warned of coarser marks since the layout */
  let warned = false;
  /** per copy its template */
  let templateOf = new Int32Array(0);
  /** per template its store indices' first and their span (a block's size: its pages' range) */
  let spanStart = new Int32Array(0);
  let spanSize = new Int32Array(0);
  /** per template its blocks let go (their offsets), reused */
  let free: number[][] = [];
  /** per copy its marker (no block), its block's offset (-1 none), the block's marked bytes */
  let uniform = new Uint8Array(0);
  let blockOf = new Int32Array(0);
  let blockMarked = new Int32Array(0);
  /**
   * per copy with a block the positions in it marked (from the block's start; a position marked again listed again,
   * cleared ones dropped on the next whole action), null where not known (a block filled whole, one or-ed in, a list
   * past a sixteenth of its size): an action on the copy whole by them alone, not its block's every byte. Why: a hover
   * cleared every copy's block whole (a template's size each: the yeast cell's 31.5 M bytes a mark)
   */
  let blockPositions: (number[] | null)[] = [];
  let words = new Uint32Array(0);
  let bytes = new Uint8Array(0);
  let used = 0;
  /**
   * the blocks' bytes changed since the last `takeChanged` (`[start, end)` pairs, in `bytes`): what the GPU's buffer
   * takes, not every block's bytes. Why: a hover uploaded every block whole (the yeast cell's copy at its levels, 36 MB
   * a mark: most of the frame)
   */
  let changedBytes: number[] = [];
  /** every byte changed (a layout: the GPU's of the one before) */
  let changedAll = false;
  const touch = (start: number, end: number) => {
    const n = changedBytes.length;
    // a position next to the last range's (a copy's particles in order) its range grown
    if (n > 0 && start >= changedBytes[n - 2] && start <= changedBytes[n - 1] + ChangedGap) {
      changedBytes[n - 1] = Math.max(changedBytes[n - 1], end);
    } else {
      changedBytes.push(start, end);
      // many: merged now (a whole action on many copies a range each)
      if (changedBytes.length > 8 * MaxChangedRanges) changedBytes = mergeChanged(changedBytes);
    }
  };

  const header = () => copyCount * 8;
  const isMarked = (c: number) => (blockOf[c] >= 0 ? blockMarked[c] > 0 : uniform[c] !== 0);
  const setWords = (c: number) => {
    const b = blockOf[c];
    if (b < 0) {
      words[2 * c] = uniform[c];
      words[2 * c + 1] = 0;
    } else {
      words[2 * c] = CopyMarkerBlock;
      words[2 * c + 1] = (header() + b - spanStart[templateOf[c]]) >>> 0;
    }
  };
  /** an action that may mark a position at 0 (the others leave it 0) */
  const isAdding = (action: MarkerAction) =>
    action === MarkerAction.Highlight || action === MarkerAction.Select || action === MarkerAction.Toggle;
  const countMarked = (from: number, to: number) => {
    let n = 0;
    for (let i = from; i < to; ++i) if (bytes[i] !== 0) n += 1;
    return n;
  };

  /** copy `c`'s block, made where it has none (filled with its marker); -1 where none fits the budget */
  function ensureBlock(c: number) {
    if (blockOf[c] >= 0) return blockOf[c];
    const t = templateOf[c];
    const size = spanSize[t];
    let b = free[t].pop();
    if (b === undefined) {
      if (used + size > budget) return -1;
      b = used;
      used += size;
      if (used > bytes.length) {
        // grown to half again within the budget (a multiple of four: the buffer's words)
        const grown = new Uint8Array(Math.ceil(Math.max(used, Math.min(budget, bytes.length * 1.5)) / 4) * 4);
        grown.set(bytes);
        bytes = grown;
      }
    }
    bytes.fill(uniform[c], b, b + size);
    // a block filled with 0 as the GPU's has it (a block let go cleared, one past the last none written)
    if (uniform[c] !== 0) touch(b, b + size);
    blockOf[c] = b;
    blockMarked[c] = uniform[c] !== 0 ? size : 0;
    blockPositions[c] = uniform[c] !== 0 ? null : [];
    setWords(c);
    return b;
  }

  /** copy `c`'s block let go where nothing in it is marked */
  function releaseIfClear(c: number) {
    const b = blockOf[c];
    if (b < 0 || blockMarked[c] > 0) return;
    free[templateOf[c]].push(b);
    blockOf[c] = -1;
    blockPositions[c] = null;
    uniform[c] = 0;
    setWords(c);
  }

  /** copy `c` marked whole for want of a block: counted, warned of once a layout */
  function noteCoarse() {
    coarseCount += 1;
    if (warned) return;
    warned = true;
    console.warn(
      `particle copies: their own markers past ${Math.round(budget / 2 ** 20)} MiB, a copy marked in part drawn marked whole`,
    );
  }

  function markCopy(c: number, action: MarkerAction) {
    if (c < 0 || c >= copyCount) return false;
    const was = isMarked(c);
    let changed: boolean;
    const b = blockOf[c];
    if (b < 0) {
      const before = uniform[c];
      applyMarkerActionAtPosition(uniform, c, action);
      changed = uniform[c] !== before;
    } else {
      const positions = blockPositions[c];
      if (positions) {
        // its marked positions alone (the others 0: an action on a whole copy marking none of them changes none, as
        // `applyMarkerAction` leaves a 0 as it is but where it adds, which a block's every byte takes: below)
        if (isAdding(action)) {
          const end = b + spanSize[templateOf[c]];
          changed = applyMarkerAction(bytes, Interval.ofBounds(b, end), action);
          if (changed) touch(b, end);
          blockMarked[c] = countMarked(b, end);
          blockPositions[c] = null;
        } else {
          changed = false;
          const kept: number[] = [];
          for (const i of positions) {
            const before = bytes[b + i];
            if (before === 0) continue;
            applyMarkerActionAtPosition(bytes, b + i, action);
            const after = bytes[b + i];
            if (after !== before) {
              changed = true;
              touch(b + i, b + i + 1);
            }
            if (after === 0) blockMarked[c] -= 1;
            else kept.push(i);
          }
          blockPositions[c] = kept;
        }
      } else {
        const end = b + spanSize[templateOf[c]];
        changed = applyMarkerAction(bytes, Interval.ofBounds(b, end), action);
        if (changed) touch(b, end);
        blockMarked[c] = countMarked(b, end);
      }
      releaseIfClear(c);
    }
    setWords(c);
    markedCount += (isMarked(c) ? 1 : 0) - (was ? 1 : 0);
    return changed;
  }

  return {
    get copyCount() {
      return copyCount;
    },
    get markedCount() {
      return markedCount;
    },
    get coarseCount() {
      return coarseCount;
    },
    get words() {
      return words;
    },
    get bytes() {
      return bytes;
    },
    get byteLength() {
      return header() + used;
    },
    takeChanged() {
      const ranges = changedAll ? undefined : mergeChanged(changedBytes);
      changedBytes = [];
      changedAll = false;
      return ranges;
    },
    layout(templates) {
      copyCount = 0;
      for (const t of templates) copyCount += t.copyCount;
      templateOf = new Int32Array(copyCount);
      spanStart = new Int32Array(templates.length);
      spanSize = new Int32Array(templates.length);
      free = templates.map(() => []);
      let c = 0;
      templates.forEach((t, k) => {
        let lo = Infinity;
        let hi = -Infinity;
        for (let i = 0; i < t.pages.length; ++i) {
          lo = Math.min(lo, t.pages[i]);
          hi = Math.max(hi, t.pages[i]);
        }
        spanStart[k] = t.pages.length ? lo * StorePageSize : 0;
        spanSize[k] = t.pages.length ? (hi - lo + 1) * StorePageSize : 0;
        templateOf.fill(k, c, c + t.copyCount);
        c += t.copyCount;
      });
      uniform = new Uint8Array(copyCount);
      blockOf = new Int32Array(copyCount).fill(-1);
      blockMarked = new Int32Array(copyCount);
      blockPositions = new Array(copyCount).fill(null);
      words = new Uint32Array(copyCount * 2);
      bytes = new Uint8Array(0);
      used = 0;
      markedCount = 0;
      coarseCount = 0;
      budget = getBudget();
      warned = false;
      // the GPU's buffer's bytes of the layout before: all of it written
      changedBytes = [];
      changedAll = true;
    },
    markCopy,
    markIndices(c, indices, action) {
      if (c < 0 || c >= copyCount || OrderedSet.size(indices) === 0) return false;
      const t = templateOf[c];
      const start = spanStart[t];
      const size = spanSize[t];
      const b = ensureBlock(c);
      if (b < 0) {
        // no block within the budget: the copy whole (where `indices` reach its template)
        let reaches = false;
        OrderedSet.forEach(indices, (index) => {
          if (index >= start && index < start + size) reaches = true;
        });
        if (!reaches) return false;
        if (isAdding(action)) noteCoarse();
        return markCopy(c, action);
      }
      const was = isMarked(c);
      let changed = false;
      OrderedSet.forEach(indices, (index) => {
        const i = index - start;
        if (i < 0 || i >= size) return;
        const before = bytes[b + i];
        applyMarkerActionAtPosition(bytes, b + i, action);
        const after = bytes[b + i];
        if (before === after) return;
        changed = true;
        touch(b + i, b + i + 1);
        blockMarked[c] += (after !== 0 ? 1 : 0) - (before !== 0 ? 1 : 0);
        const positions = blockPositions[c];
        if (before === 0 && positions) {
          if (positions.length * 16 < size) positions.push(i);
          else blockPositions[c] = null;
        }
      });
      releaseIfClear(c);
      markedCount += (isMarked(c) ? 1 : 0) - (was ? 1 : 0);
      return changed;
    },
    markAll(action) {
      let changed = false;
      for (let c = 0; c < copyCount; ++c) if (markCopy(c, action)) changed = true;
      return changed;
    },
    include(c, values) {
      if (c < 0 || c >= copyCount) return false;
      const t = templateOf[c];
      const start = spanStart[t];
      const size = spanSize[t];
      let any = false;
      for (let i = 0; i < size && !any; ++i) any = values[start + i] !== 0;
      if (!any) return false;
      const was = isMarked(c);
      const b = ensureBlock(c);
      if (b < 0) {
        // no block within the budget: the copy whole, its marker every value's
        let v = 0;
        for (let i = 0; i < size; ++i) v |= values[start + i];
        noteCoarse();
        const before = uniform[c];
        uniform[c] |= v;
        setWords(c);
        markedCount += (isMarked(c) ? 1 : 0) - (was ? 1 : 0);
        return uniform[c] !== before;
      }
      blockPositions[c] = null;
      touch(b, b + size);
      for (let i = 0; i < size; ++i) {
        const v = values[start + i];
        if (v === 0) continue;
        const before = bytes[b + i];
        bytes[b + i] = before | v;
        if (before === 0) blockMarked[c] += 1;
      }
      markedCount += (isMarked(c) ? 1 : 0) - (was ? 1 : 0);
      return true;
    },
    isMarked: (c) => c >= 0 && c < copyCount && isMarked(c),
    getMarkedPositions: (c) => (blockOf[c] >= 0 ? blockPositions[c] : undefined),
    get(c, index) {
      const b = blockOf[c];
      if (b < 0) return uniform[c];
      const i = index - spanStart[templateOf[c]];
      return i >= 0 && i < spanSize[templateOf[c]] ? bytes[b + i] : 0;
    },
  };
}
