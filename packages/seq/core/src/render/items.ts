/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn, type Column } from '@molstar/core/repr/tree';
import type { DataEntry, RangesDensity, RangesList } from '@molstar/seqx-core/data';
import { rangesLocation } from '@molstar/seqx-core/query';
import { type LayerSource, toSection } from './layer-data.js';
import { addToStrip, paintStrip, resetStrip, type StripBuffers } from './ranges.js';
import type { SeqTokens } from './tokens.js';
import type { AxisTransform } from './axis.js';
import type { Box, Hit } from './layer-type.js';

// The items of a ranges layer's tiles for the layers that draw them their own way (`secondary`, `variants`; _spec/README.md
// "Render", Layer types, Data): each row an item with its key, extent and the value of one field; a key and start seen
// in an earlier tile is that range again (over the tile's edge), a key repeated with another start is another piece of
// the item.

export interface ItemRows {
  /** The rows in view (`n` of them; the arrays are kept longer and written over frame to frame). */
  n: number;
  keys: (string | number)[];
  starts: number[];
  ends: number[];
  /** The value of the field asked for, as text. */
  values: (string | undefined)[];
}

export function newItemRows(): ItemRows {
  return { n: 0, keys: [], starts: [], ends: [], values: [] };
}

/** A table's field by name. */
function fieldOf(fields: RangesList['fields'], name: string): Column | undefined {
  if (!fields) return undefined;
  for (let k = 0; k < fields.length; k++) if (fields[k].name === name) return fields[k].column;
  return undefined;
}

function textAt(c: Column | undefined, i: number): string | undefined {
  if (!c) return undefined;
  const v = c instanceof CategoryColumn ? c.get(i) : c[i];
  return v === null || v === undefined ? undefined : String(v);
}

export interface DensityTile {
  readonly table: RangesDensity;
  readonly level: number;
}

/**
 * Reads the list tables of `entries` into `rows` (cleared first), each row's value of the first of `fields` the table
 * has; the density tables are returned.
 */
export function collectItems(
  rows: ItemRows,
  entries: readonly Extract<DataEntry, { status: 'loaded' }>[],
  fields: readonly string[],
): DensityTile[] {
  rows.n = 0;
  const densities: DensityTile[] = [];
  const earlier = new Map<string | number, number[]>();
  for (const e of entries) {
    const t = e.tables.ranges;
    if (!t) continue;
    if (t.form === 'density') {
      densities.push({ table: t, level: e.region.level });
      continue;
    }
    const from = rows.n;
    let col: Column | undefined;
    for (let f = 0; f < fields.length && !col; f++) col = fieldOf(t.fields, fields[f]);
    for (let i = 0; i < t.start.length; i++) {
      const id = t.id?.get(i);
      const key: string | number =
        id !== null && id !== undefined ? id : t.row ? t.row[i] : `${e.region.level}/${e.region.tile}:${i}`;
      if (earlier.get(key)?.includes(t.start[i])) continue;
      const j = rows.n++;
      rows.keys[j] = key;
      rows.starts[j] = t.start[i];
      rows.ends[j] = t.end[i];
      rows.values[j] = textAt(col, i);
    }
    for (let j = from; j < rows.n; j++) {
      const starts = earlier.get(rows.keys[j]);
      if (starts) starts.push(rows.starts[j]);
      else earlier.set(rows.keys[j], [rows.starts[j]]);
    }
  }
  return densities;
}

/**
 * An item's hit: its element with its extent on the data's system and its hull on the section's (every piece carried
 * along the layer's path, as the ranges layer's), the position under the pointer.
 */
export function itemHit(
  base: Pick<Hit, 'layer' | 'section' | 'system' | 'position'>,
  owner: string,
  dataSystem: string,
  src: LayerSource | undefined,
  key: string | number,
  s: number,
  e: number,
): Hit {
  let lo = s;
  let hi = e;
  if (src?.path) {
    const pieces = toSection(src, s, e, []);
    lo = Number.POSITIVE_INFINITY;
    hi = Number.NEGATIVE_INFINITY;
    for (let k = 0; k + 1 < pieces.length; k += 2) {
      lo = Math.min(lo, pieces[k]);
      hi = Math.max(hi, pieces[k + 1]);
    }
  }
  return {
    ...base,
    element: { owner, key },
    extent: rangesLocation(dataSystem, [s, e]),
    ...(lo < hi ? { hull: [lo, hi] as const } : {}),
  };
}

/**
 * Which rows hold the ends of their item, in reading order (rows sharing a key are the pieces of one item): `last` the
 * row holding its end, `first` the one holding its start. Across the origin of a cyclic system the row that starts at it
 * comes last (the longest, if several) and the one that ends there first; otherwise the row ending last and the row
 * starting first. `first` and `last` are written for rows `0..n`.
 */
export function markEnds(
  keys: ArrayLike<string | number>,
  starts: ArrayLike<number>,
  ends: ArrayLike<number>,
  n: number,
  cyclic: boolean,
  length: number,
  first: Uint8Array,
  last: Uint8Array,
) {
  const groups = new Map<string | number, number[]>();
  for (let i = 0; i < n; i++) {
    first[i] = 0;
    last[i] = 0;
    const g = groups.get(keys[i]);
    if (g) g.push(i);
    else groups.set(keys[i], [i]);
  }
  for (const g of groups.values()) {
    const across = cyclic && g.length > 1 && g.some((i) => ends[i] === length) && g.some((i) => starts[i] === 0);
    let c = across ? -1 : g[0];
    let f = across ? -1 : g[0];
    for (const i of g) {
      if (across) {
        if (starts[i] === 0 && (c < 0 || ends[i] > ends[c])) c = i;
        if (ends[i] === length && (f < 0 || starts[i] < starts[f])) f = i;
      } else {
        if (ends[i] > ends[c]) c = i;
        if (starts[i] < starts[f]) f = i;
      }
    }
    last[c] = 1;
    first[f] = 1;
  }
}

const imageEnd: number[] = [];
const imageStart: number[] = [];

/**
 * Where the first and last positions of a row `[s, e)` land on the section, and whether the row runs against it (its
 * start image above its end image): what an arrowhead at the row's end or start needs. On the section's own system
 * the positions themselves. No number made per row: read the fields, never a call with them.
 */
export class RowImages {
  /** The last position's image `[ea, eb)` and the first's `[sa, sb)`; NaN when a position maps nowhere. */
  ea = Number.NaN;
  eb = Number.NaN;
  sa = Number.NaN;
  sb = Number.NaN;
  flip = false;

  set(src: LayerSource | undefined, s: number, e: number) {
    this.flip = false;
    if (!src?.path) {
      this.ea = e - 1;
      this.eb = e;
      this.sa = s;
      this.sb = s + 1;
      return;
    }
    toSection(src, e - 1, e, imageEnd);
    toSection(src, s, s + 1, imageStart);
    if (imageEnd.length < 2 || imageStart.length < 2) {
      this.ea = this.eb = this.sa = this.sb = Number.NaN;
      return;
    }
    this.ea = imageEnd[0];
    this.eb = imageEnd[imageEnd.length - 1];
    this.sa = imageStart[0];
    this.sb = imageStart[imageStart.length - 1];
    this.flip = e - s > 1 && this.sa > this.ea;
  }

  /**
   * The side an arrowhead stands on at the end (`atEnd`) or the start of the row, for the piece `[pa, pb)` of the
   * section it is drawn on, as it is drawn: 0 none (the end is not in this piece, or the view clips it), 1 the tip on
   * the right, 2 on the left. The end's side is the screen's right where the residues run to the right.
   */
  tip(axis: AxisTransform, box: Box, pa: number, pb: number, atEnd: boolean): number {
    const lo = atEnd ? this.ea : this.sa;
    const hi = atEnd ? this.eb : this.sb;
    if (!(lo >= pa && hi <= pb)) return 0;
    const at = atEnd === this.flip ? lo : hi;
    const x = axis.x(at);
    if (!(x >= box.x - 0.5 && x <= box.x + box.width + 0.5)) return 0;
    return ((axis.dir === 1) !== this.flip) === atEnd ? 1 : 2;
  }
}

/**
 * The count strip of the items and the density bins of a read, one a device column, in the density tokens: its image
 * (the layers' strip at the `bars` and `density` steps and from a source's bins).
 */
export function countStrip(
  strip: StripBuffers,
  src: LayerSource,
  axis: AxisTransform,
  rows: ItemRows,
  densities: readonly DensityTile[],
  tokens: SeqTokens,
): StripBuffers['image'] {
  resetStrip(strip, axis.columns);
  strip.v = 1;
  for (let i = 0; i < rows.n; i++)
    addToStrip(strip, src, axis, rows.starts[i], Math.max(rows.ends[i], rows.starts[i] + 1));
  for (const { table, level } of densities)
    for (let i = 0; i < table.bin.length; i++) {
      strip.v = table.count[i];
      addToStrip(strip, src, axis, table.bin[i] * level, (table.bin[i] + 1) * level);
    }
  return paintStrip(strip, axis.columns, tokens);
}
