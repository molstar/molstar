/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `color` renderer (seqspec:viewer#extension-points): a cell per position by the color theme (wrapped:
// backgrounds), reduced by the theme's kind below a device px (seqspec:sources#levels-and-the-zoom-gate): a categorical
// theme by its classes, keyed by row where it colors per row or coarser; a continuous one by `palette(mean)` over the
// dense values' prepared pyramid or a row field's values. A position no row covers and with no value is none of the
// layer's.
import {
  type CellSource,
  cellBand,
  type DrawView,
  drawCells,
  Granularity,
  kernelEntry,
  NO_COLOR,
  OPAQUE,
  type PerPosition,
  positionKey,
  RectBatch,
  rowValues,
  ThemeLocation,
  visibleRange,
} from '../draw/api.js';

const rects = new RectBatch();
const range = new Int32Array(2);

const src: CellSource & { v: DrawView | null; s: PerPosition | null } = {
  v: null,
  s: null,
  color(b, i) {
    const v = src.v as DrawView;
    const s = src.s as PerPosition;
    if (!v.color || !v.axis) return NO_COLOR;
    const j = i - v.axis.base[b];
    const rows = s.row[b];
    const values = s.value[b];
    const r = rows ? rows[j] : -1;
    const value = values ? values[j] : Number.NaN;
    if (r < 0 && Number.isNaN(value)) return NO_COLOR;
    const ann = s.annotation[b];
    const loc = ThemeLocation.at(v.loc, v.axis, v.layer, b, i);
    loc.value = value;
    loc.annotation = ann ? ann[j] : s.only[b];
    loc.row = r;
    return v.color.value(loc);
  },
  alpha() {
    const t = (src.v as DrawView).color;
    return t?.alpha ? t.alpha((src.v as DrawView).loc) : OPAQUE;
  },
  key(b, i) {
    const v = src.v as DrawView;
    const g = v.color ? v.color.granularity : -1;
    const byRow =
      g === Granularity.unit || g === Granularity.row || g === Granularity.item || g === Granularity.uniform;
    return positionKey(src.s as PerPosition, b, i - (v.axis ? v.axis.base[b] : 0), byRow);
  },
  // A theme over the values reads them (and their prepared pyramid); one over a row field, each row's value.
  meanValues(b) {
    const t = (src.v as DrawView).color;
    const s = src.s as PerPosition;
    if (t?.granularity === Granularity.value || !t?.rowValue) return s.value[b];
    return s.row[b] ? rowValues(s, b, t.rowValue) : null;
  },
  meanPyramid(b) {
    const t = (src.v as DrawView).color;
    return t?.granularity === Granularity.value || !t?.rowValue ? (src.s as PerPosition).pyramid[b] : null;
  },
};

function draw(v: DrawView, s: PerPosition) {
  rects.reset();
  ThemeLocation.reset(v.loc);
  const pad = v.height > 8 ? 2 : 0;
  cellBand[0] = Math.round((v.top + pad) * v.dpr);
  cellBand[1] = Math.round((v.top + v.height - pad) * v.dpr) - cellBand[0];
  cellBand[2] = 0;
  src.v = v;
  src.s = s;
  try {
    for (let b = 0; b < s.value.length; b++) {
      visibleRange(v, b, range);
      if (s.value[b] || s.row[b]) drawCells(v, b, range[0], range[1], src, s, rects);
    }
    if (rects.n > 0) v.backend.rects(rects);
  } finally {
    src.v = src.s = null;
  }
}

export const ColorRenderer = kernelEntry<PerPosition>('color', 14, {
  geometries: ['linear', 'wrapped'],
  stable: true,
  draw,
  tiny: draw,
});
