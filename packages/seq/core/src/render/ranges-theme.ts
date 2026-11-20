/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { CategoryColumn } from '@molstar/core/repr/tree';
import type { ItemRows } from './items.js';
import type { Channel } from '@molstar/seqx-core/theme/index';
import type { LayoutContext } from './layer-type.js';
import { rangesColorFallback } from './layer-color.js';
import { withAlpha } from './tokens.js';
import type { RangesLay, RangesPrep } from './ranges-lay.js';

// A ranges layer's channels (_spec/README.md "Themes"): each item's fill, size and painted value through its table's
// row, a colourless layer's categorical default, a dictionary-coded column's colors once a table.

/**
 * The items' channels (_spec/README.md "Themes"): each item's fill, size and painted value, read through its table's
 * row. The channels are compiled once a layout (kept while the tokens and held domains agree), by the first table's
 * declared column types; a field theme over a dictionary-coded column is a color per code, once a table.
 */
export function themeItems(lay: RangesLay, prep: RangesPrep, ctx: LayoutContext) {
  const rows = lay.rows;
  const datum = lay.datum;
  const missing = ctx.tokens.colors['seq-missing'];
  datum.bind(rows.n ? rows.fields[0] : undefined);
  const extra = datum.env;
  // a layer with no color: the next categorical color by its order in the view (_spec/README.md "Themes")
  const color = prep.channels.get('color', rangesColorFallback(ctx.colorIndex), ctx, extra);
  const size = prep.channels.get('size', 1, ctx, extra);
  const opacity = prep.channels.get('opacity', 1, ctx, extra);
  lay.color = color;
  lay.size = size;
  lay.extents.reset([color, size, opacity]);
  const painted = color.legend?.kind === 'palette' ? color.legend.field : undefined;
  const cconst = color.constant;
  const sconst = size.constant;
  const oconst = opacity.constant;
  let bound: ItemRows['fields'][number] | null = null;
  let codes: ArrayLike<number> | undefined;
  let dict: readonly string[] | undefined;
  for (let j = 0; j < rows.n; j++) {
    const fields = rows.fields[j];
    if (fields !== bound) {
      bound = fields;
      datum.bind(fields);
      const coded = cconst === undefined && color.field !== undefined ? datum.column(color.field) : undefined;
      codes = coded instanceof CategoryColumn ? coded.codes : undefined;
      dict = coded instanceof CategoryColumn ? coded.dict : undefined;
      if (dict) lay.byCode = codeFills(color, dict, missing, lay.byCode);
    }
    const i = rows.row[j];
    datum.row = i;
    datum.start = rows.starts[j];
    datum.end = rows.ends[j];
    let fill: number;
    if (codes) {
      const code = codes[i];
      fill = code < 0 ? missing : lay.byCode[code];
    } else {
      const c = cconst ?? color.at(datum);
      fill = typeof c === 'number' ? c : missing;
    }
    const o = oconst ?? opacity.at(datum);
    if (typeof o === 'number' && o < 1) fill = withAlpha(fill, ((fill & 0xff) / 255) * Math.max(0, o));
    const z = sconst ?? size.at(datum);
    lay.sizes[j] = typeof z === 'number' ? Math.max(0, Math.min(1, z)) : 1;
    if (painted !== undefined) {
      if (codes) lay.painted[j] = codes[i] < 0 ? undefined : dict![codes[i]];
      else {
        const v = datum.field(painted);
        lay.painted[j] = typeof v === 'string' ? v : undefined;
      }
    }
    if (lay.extents.any) lay.extents.add(datum);
    lay.fills[j] = fill;
  }
  lay.extents.flush(ctx.theme);
}

/** The colors of a dictionary's values through a channel that reads one field: index a code. */
function codeFills(channel: Channel, dict: readonly string[], missing: number, out: Uint32Array): Uint32Array {
  if (out.length < dict.length) out = new Uint32Array(Math.max(dict.length, 2 * out.length));
  let v: string | undefined;
  const d = { value: Number.NaN, letter: '', start: 0, end: 0, field: () => v };
  for (let k = 0; k < dict.length; k++) {
    v = dict[k];
    const c = channel.at(d);
    out[k] = typeof c === 'number' ? c : missing;
  }
  return out;
}
