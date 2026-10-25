/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author Alexander Rose <alexander.rose@weirdbyte.de>
 */

import { ValueCell } from '@molstar/core/util/value-cell';
import type { DataArray, TextureImage } from './util.js';

/**
 * The element ranges the updates of a data array (`DataSpec`) changed, so that consumers upload those instead of the whole image
 * (`shared-draw-resources.md`): a short log on the image, an entry per update by the cell version it produced, its
 * ranges as `[start, end)` element pairs. Only `TextureChanges.update` writes entries, together with the update; any
 * other update of the cell has none, and consumers upload all then. Changes of the image that do not update the cell
 * (`defer`) are pending, part of the next update's entry
 */
type ChangeLog = {
  readonly versions: number[];
  readonly ranges: (number[] | undefined)[];
  /** the ranges changed since the last update, not uploaded yet (`defer`); undefined: none, `null`: all */
  pending?: number[] | null;
};

const logs = new WeakMap<DataArray<any>, ChangeLog>();

/** Updates kept per image: consumers further behind upload all */
const MaxEntries = 8;
/** Ranges per update kept: more upload all */
const MaxEntryRanges = 64;
/** Ranges closer than this (in elements) are merged: uploads of what lies between cost less than another write */
const MergeGap = 256;
/** Ranges handed out at most, and the fill (of all elements) beyond which consumers upload all */
const MaxRanges = 8;
const MaxFill = 0.5;

/** The rows each in-place write of an image changed, by the image version it produced (`TextureChanges.writeRows`) */
type RowLog = { readonly versions: number[]; readonly rows: [number, number][] };
const rowLogs = new WeakMap<TextureImage<any>, RowLog>();
/** Writes kept per image (a label's new glyphs, a write each, between two uploads): consumers further behind upload all */
const MaxRowEntries = 64;

export namespace TextureChanges {
  /**
   * Updates `cell` (its image changed in place) with the element ranges that changed (`[start, end)` pairs), logged
   * for the version the update produces; `ranges` undefined: the whole image changed
   */
  export function update(cell: ValueCell<DataArray<any>>, ranges: readonly number[] | undefined) {
    const image = cell.ref.value;
    let log = logs.get(image);
    if (!log) logs.set(image, (log = { versions: [], ranges: [] }));
    const { pending } = log;
    log.pending = undefined;
    const all = !ranges || pending === null;
    const entry = all ? undefined : pending ? [...pending, ...ranges] : ranges.slice();
    log.versions.push(cell.ref.version + 1);
    log.ranges.push(entry && entry.length <= 2 * MaxEntryRanges ? entry : undefined);
    if (log.versions.length > MaxEntries) {
      log.versions.shift();
      log.ranges.shift();
    }
    ValueCell.update(cell, image);
  }

  /**
   * Notes that the element ranges of the image of `cell` changed without an update of the cell (e.g. what a uniform
   * stands for meanwhile, as the marker status): the next `update` uploads them too; `ranges` undefined: all of it.
   * Without, consumers holding the image would keep what they uploaded before for these ranges
   */
  export function defer(cell: ValueCell<DataArray<any>>, ranges: readonly number[] | undefined) {
    const image = cell.ref.value;
    let log = logs.get(image);
    if (!log) logs.set(image, (log = { versions: [], ranges: [] }));
    if (!ranges || log.pending === null) log.pending = null;
    else log.pending = log.pending ? [...log.pending, ...ranges] : ranges.slice();
    if (log.pending && log.pending.length > 2 * MaxEntryRanges) log.pending = null;
  }

  /**
   * The element ranges of `image` changed from version `from` to `to` of its cell, merged and sorted: undefined where
   * an update in between is not logged, changed all, or the ranges are many or fill much of `elementCount`
   */
  export function get(image: DataArray<any>, from: number, to: number, elementCount = image.count) {
    const log = logs.get(image);
    if (!log || to <= from || to - from > log.versions.length) return undefined;
    const pairs: number[] = [];
    for (let v = from + 1; v <= to; ++v) {
      const i = log.versions.lastIndexOf(v);
      const ranges = i < 0 ? undefined : log.ranges[i];
      if (!ranges) return undefined;
      for (let k = 0; k < ranges.length; k += 2) pairs.push(ranges[k], ranges[k + 1]);
    }
    return mergeRanges(pairs, elementCount);
  }

  /**
   * The rows (`[start, end)` pairs, merged) of `image` changed from version `from` to `to` of its cell, laid out in rows
   * of `width` elements (a data texture): undefined where `get` is
   */
  export function getRows(image: DataArray<any>, from: number, to: number, width: number) {
    const ranges = get(image, from, to);
    if (!ranges) return undefined;
    const rows: number[] = [];
    for (let k = 0; k < ranges.length; k += 2) {
      const start = Math.floor(ranges[k] / width);
      const end = Math.ceil(ranges[k + 1] / width);
      const last = rows.length - 1;
      if (last > 0 && start <= rows[last]) rows[last] = Math.max(rows[last], end);
      else rows.push(start, end);
    }
    return rows;
  }

  /**
   * Notes that `image` changed in place in its rows `[start, end)` (a glyph appended to a font atlas): its `version`
   * advances, logged with the rows, so that a texture holding an earlier version uploads those rows (`getImageRows`,
   * `Texture.loadRows`) instead of all of it
   */
  export function writeRows(image: TextureImage<any>, start: number, end: number) {
    const version = (image.version ?? 0) + 1;
    let log = rowLogs.get(image);
    if (!log) rowLogs.set(image, (log = { versions: [], rows: [] }));
    log.versions.push(version);
    log.rows.push([start, end]);
    if (log.versions.length > MaxRowEntries) {
      log.versions.shift();
      log.rows.shift();
    }
    (image as { version?: number }).version = version;
  }

  /**
   * The rows (`[start, end)` pairs, sorted and merged) of `image` written from its version `from` to `to`
   * (`writeRows`): undefined where a version in between is not logged, or they cover much of the image (`MaxFill`)
   */
  export function getImageRows(image: TextureImage<any>, from: number, to: number) {
    const log = rowLogs.get(image);
    if (!log || to <= from || to - from > log.versions.length) return undefined;
    const written: [number, number][] = [];
    for (let v = from + 1; v <= to; ++v) {
      const i = log.versions.lastIndexOf(v);
      if (i < 0) return undefined;
      written.push(log.rows[i]);
    }
    written.sort((a, b) => a[0] - b[0]);
    const rows: number[] = [];
    let covered = 0;
    for (const [s, e] of written) {
      const start = Math.max(0, s);
      const end = Math.min(image.height, e);
      if (end <= start) continue;
      const last = rows.length - 1;
      if (last > 0 && start <= rows[last]) rows[last] = Math.max(rows[last], end);
      else rows.push(start, end);
    }
    for (let k = 0; k < rows.length; k += 2) covered += rows[k + 1] - rows[k];
    return covered > MaxFill * image.height ? undefined : rows;
  }

  /**
   * `pairs` sorted and merged where closer than `MergeGap`, clamped to `elementCount`: undefined past `MaxRanges` or a
   * fill of `MaxFill` (one upload of all costs less than many scattered ones)
   */
  export function mergeRanges(pairs: readonly number[], elementCount: number) {
    const order: number[] = [];
    for (let k = 0; k < pairs.length; k += 2) order.push(k);
    order.sort((a, b) => pairs[a] - pairs[b]);
    const merged: number[] = [];
    let covered = 0;
    for (const k of order) {
      const start = Math.max(0, pairs[k]);
      const end = Math.min(elementCount, pairs[k + 1]);
      if (end <= start) continue;
      const last = merged.length - 1;
      if (last > 0 && start - merged[last] < MergeGap) {
        if (end > merged[last]) merged[last] = end;
      } else {
        merged.push(start, end);
      }
    }
    if (merged.length > 2 * MaxRanges) return undefined;
    for (let k = 0; k < merged.length; k += 2) covered += merged[k + 1] - merged[k];
    if (covered > MaxFill * elementCount) return undefined;
    return merged;
  }
}
