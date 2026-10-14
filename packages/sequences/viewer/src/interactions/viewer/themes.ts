/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `interactions:kind` theme's implementation (seqspec:viewer#extension-points; seqspec:libraries#interactions): an
// interaction's kind through the vocabulary's colors and shapes, the same ones the format half writes out as kernel
// fallbacks. One provider serves the color and the glyph channel, as its one format entry does. A location's kind is
// its `category` where a renderer sets it (the fingerprint, per cell and kind: a kind index), else its pair row's
// category. Tables are built in `create`; `value` reads them, allocating nothing.

import { Color } from '@molstar/seq-core/format';
import {
  dataIndex,
  GlyphShape,
  Granularity,
  NO_COLOR,
  NO_GLYPH,
  type Theme,
  type ThemeContext,
  type ThemeLocation,
  type ThemeProvider,
} from '@molstar/seq-viewer/draw';
import { fingerprintOf } from '../registry.js';
import { KindStyle, Kinds, kindByName } from '../vocabulary.js';

const COLORS = Int32Array.from(Kinds, (k) => Color.decode(KindStyle[k][1]) ?? NO_COLOR);
const SHAPES = Int32Array.from(Kinds, (k) => GlyphShape[KindStyle[k][2] as keyof typeof GlyphShape] ?? NO_GLYPH);

/** Per data annotation, per kept pair row: its kind (-1 elsewhere); and the kinds present (fingerprints' cells too). */
function kindsOf(tc: ThemeContext) {
  let present = 0;
  const rows = tc.annotations.map((a, k) => {
    const kept = tc.layer.kept[k];
    const fp = a ? fingerprintOf(a) : undefined;
    for (let c = 0; fp && c < fp.cells.length; c++) present |= fp.cells[c];
    const out = new Int8Array(a?.count ?? 0).fill(-1);
    for (let r = 0; a?.pairs && r < a.count; r++) {
      if (kept && !kept[r]) continue;
      out[r] = kindByName(a.category && a.category[r] >= 0 ? tc.strings.get(a.category[r]) : null);
      present |= 1 << out[r];
    }
    return out;
  });
  /** The location's kind; -1: none. */
  const kindAt = (loc: ThemeLocation) => {
    if (loc.category >= 0) return loc.category < Kinds.length ? loc.category : -1;
    const k = loc.annotation >= 0 && loc.item >= 0 ? dataIndex(tc, loc.annotation) : -1;
    return k >= 0 && loc.item < rows[k].length ? rows[k][loc.item] : -1;
  };
  return { present, kindAt };
}

export const KindTheme: ThemeProvider = {
  name: 'interactions:kind',
  label: 'Interaction kind',
  category: 'interactions',
  kinds: ['color', 'glyph'],
  granularity: Granularity.item,
  isApplicable: () => true,
  create(tc, _params, channel): Theme {
    const { present, kindAt } = kindsOf(tc);
    const color = channel.kind === 'color';
    const table = color ? COLORS : SHAPES;
    const none = color ? NO_COLOR : NO_GLYPH;
    return {
      kind: channel.kind,
      granularity: Granularity.item,
      reduction: 'classes',
      value: (loc) => {
        const k = kindAt(loc);
        return k >= 0 ? table[k] : none;
      },
      contextHash: tc.stamp(),
      /** The kinds present, in vocabulary order, each with its color and shape. */
      legend: color
        ? (out) => {
            for (let k = 0; k < Kinds.length; k++)
              if (present & (1 << k)) out.item(KindStyle[Kinds[k]][0], COLORS[k], SHAPES[k]);
          }
        : undefined,
    };
  },
};
