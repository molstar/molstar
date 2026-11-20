/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { grown, NoColor } from './draw-lists.js';
import { ImageTile } from './images.js';
import { rowsLevel, rowsPending } from './rows-levels.js';
import type { LayoutContext } from './layer-type.js';
import { NoClass, Pyramid, pyramidLevelFor } from './pyramid.js';
import { consensus, type HeldCells, upper } from './rows-model.js';
import { inkOn, type SeqTokens } from './tokens.js';
import type { RowsLay, RowsPrep } from './rows-lay.js';

// The rows layer's cells in view (_spec/README.md "Render", The rows layer): rects (letters at the letters step), a
// folded clade's consensus, or below 3 device px the shared reduction's image tiles.

/**
 * The cells in view: rects (letters at the letters step), or image tiles from the reduction below 3 device px; the
 * pinned rows after the scrolled ones (`pinnedFrom`), drawn on their own over a band that does not scroll.
 */
export function layoutCells(lay: RowsLay, prep: RowsPrep, ctx: LayoutContext) {
  const { axis, config, tokens } = ctx;
  const dpr = axis.dpr;
  const colPx = axis.devicePx;
  const rowPx = lay.rowScale * dpr;
  const vis = axis.visible([]);
  if (!vis.length) return;
  const c0 = Math.max(0, vis[0]);
  const c1 = Math.min(lay.columns, vis[vis.length - 1]);
  lay.viewC0 = c0;
  lay.viewC1 = c1;
  const missing = tokens.colors['seq-missing'];
  const held = lay.held;
  const C = lay.columns;
  if ((colPx < config.density.cells || rowPx < config.density.cells) && held) layoutImages(lay, ctx, c0, c1);
  else {
    const lettersStep = axis.scale >= config.density.letters && lay.rowScale >= config.density.letters;
    lay.letters = lettersStep;
    if (lettersStep) {
      let size = 0;
      for (const s of config.letterSizes)
        if (s * config.monoAdvance <= axis.scale - 2 && s <= lay.rowScale * 0.75) size = s;
      lay.fontSize = size || (config.letterSizes[0] ?? 9);
      lay.font = `${lay.fontSize}px mono`;
    }
    for (let i = lay.first; i < lay.last; i++) addRow(lay, i, c0, c1, held, C, missing, tokens, prep);
  }
  lay.pinnedFrom = lay.n;
  for (let i = 0; i < lay.display.pinned; i++) addRow(lay, i, c0, c1, held, C, missing, tokens, prep);
}

function addRow(
  lay: RowsLay,
  i: number,
  c0: number,
  c1: number,
  held: HeldCells | undefined,
  C: number,
  missing: number,
  tokens: SeqTokens,
  prep: RowsPrep,
) {
  const axis = lay.axis;
  const dpr = axis.dpr;
  const e = lay.display.rowAt(i);
  const y = lay.rowY(i);
  const h = lay.rowScale;
  let codes: Uint8Array | undefined;
  let base = 0;
  let summary: Uint8Array | undefined;
  if (e < 0) summary = summaryOf(lay, i);
  else if (held) {
    codes = held.codes;
    base = e * C;
  } else codes = lay.lazyCodes.get(e);
  if (!codes && !summary) return;
  let run = -1;
  for (let j = c0; j < c1; j++) {
    let code: number;
    let alpha = 255;
    if (summary) {
      code = summary[3 * j];
      if (prep.summary === 'bar') {
        alpha = summary[3 * j + 2];
        code = alpha ? 65 : 45;
      } else alpha = summary[3 * j + 1];
    } else code = codes![base + j];
    let fill = lay.colorOf(code, j, missing);
    if (fill === NoColor) {
      run = -1;
      continue;
    }
    if (alpha < 255) fill = ((fill & 0xffffff00) | Math.round(((fill & 0xff) * alpha) / 255)) >>> 0;
    if (!axis.cell(j)) continue;
    const x0 = Math.round(axis.sx0 * dpr) / dpr;
    const x1 = Math.round(axis.sx1 * dpr) / dpr;
    // a run of one fill on a row at the cells step is one rect (bounded by the screen)
    if (!lay.letters && run >= 0 && lay.fill[run] === fill && lay.x1[run] >= x0 - 0.01) {
      lay.x1[run] = x1;
      continue;
    }
    const k = lay.n++;
    if (k >= lay.x0.length) {
      lay.x0 = grown(lay.x0, k + 1);
      lay.x1 = grown(lay.x1, k + 1);
      lay.y = grown(lay.y, k + 1);
      lay.h = grown(lay.h, k + 1);
      lay.fill = grown(lay.fill, k + 1);
      lay.code = grown(lay.code, k + 1);
      lay.ink = grown(lay.ink, k + 1);
    }
    lay.x0[k] = x0;
    lay.x1[k] = x1;
    lay.y[k] = y;
    lay.h[k] = h;
    lay.fill[k] = fill;
    lay.code[k] = summary && prep.summary === 'bar' ? 0 : code;
    lay.ink[k] = lay.letters ? inkOn(fill, tokens) : NoColor;
    run = k;
  }
}

/** A folded clade's consensus at display row `i`, computed exactly from its cells when first drawn and kept. */
function summaryOf(lay: RowsLay, i: number): Uint8Array | undefined {
  const held = lay.held;
  const members = lay.display.members.get(i);
  if (!held || !members) return undefined;
  if (lay.display.version !== lay.summariesOf) {
    lay.summaries.clear();
    lay.summariesOf = lay.display.version;
  }
  let s = lay.summaries.get(i);
  if (!s) lay.summaries.set(i, (s = consensus(held, members, new Uint8Array(3 * held.columns))));
  return s;
}

/**
 * Below 3 device px a cell: the shared reduction over the scrolled display rows' theme classes (a class a distinct
 * fill), at the level pair whose sample is at least a device pixel, drawn as image tiles (256 × 256 samples) kept in
 * the shared image cache by (layer, level, tile, base, theme), offset from display row `pinned` (the pinned rows are
 * drawn on their own).
 */
function layoutImages(lay: RowsLay, ctx: LayoutContext, c0: number, c1: number) {
  const { axis, tokens } = ctx;
  const held = lay.held!;
  const display = lay.display;
  const C = held.columns;
  const P = display.pinned;
  const pyr = (lay.pyramid ??= new Pyramid({ factor: 4, pairing: 'independent', reducer: 'majority' }));
  const baseKey = `${held.version}|${display.version}|${tokens.version}|${lay.channelVersion}`;
  const R = display.count - P;
  if (lay.pyramidOf !== baseKey) {
    lay.pyramidOf = baseKey;
    const need = R * C;
    if (!lay.baseBytes || lay.baseBytes.length < need) lay.baseBytes = new Uint8Array(need);
    const base = lay.baseBytes;
    // one class a distinct fill (at most 254): the theme's classes, so a reduced sample reads as the theme; a table by
    // (column, letter) made once, then one lookup a cell
    const classOf = new Map<number, number>();
    const missing = tokens.colors['seq-missing'];
    const table = lay.classTable && lay.classTable.length >= C * 128 ? lay.classTable : new Uint8Array(C * 128);
    lay.classTable = table;
    table.fill(254, 0, C * 128);
    const classAt = (code: number, j: number): number => {
      const u = upper(code);
      if (u > 127) return NoClass;
      const k = j * 128 + u;
      let c = table[k];
      if (c !== 254) return c;
      const fill = lay.colorOf(code, j, missing);
      if (fill === NoColor) c = NoClass;
      else {
        let id = classOf.get(fill);
        if (id === undefined) {
          id = classOf.size < 253 ? classOf.size : 252;
          if (classOf.size < 253) classOf.set(fill, id);
          lay.classColors[id] = fill;
        }
        c = id;
      }
      table[k] = c;
      return c;
    };
    for (let i = 0; i < R; i++) {
      const e = display.rowAt(P + i);
      const sum = e < 0 ? summaryOf(lay, P + i) : undefined;
      if (sum) {
        for (let j = 0; j < C; j++) base[i * C + j] = classAt(sum[3 * j], j);
        continue;
      }
      const row = e * C;
      for (let j = 0; j < C; j++) {
        const code = held.codes[row + j];
        const c = table[j * 128 + (code & 127)];
        base[i * C + j] = c !== 254 && code < 97 ? c : classAt(code, j);
      }
    }
    pyr.setClasses(R, C, base.subarray(0, need));
  }
  const dpr = axis.dpr;
  const rowLevel = pyramidLevelFor(1 / (lay.rowScale * dpr), 4);
  const colLevel = pyramidLevelFor(1 / axis.devicePx, 4);
  const level = R > 0 ? rowsLevel(pyr, rowLevel, colLevel, ctx.settled, ctx.config.rows.buildBudget) : undefined;
  if (!level) return void (R > 0 && rowsPending(lay));
  const fr = 4 ** level.rowLevel;
  const fc = 4 ** level.columnLevel;
  const s0 = Math.floor((lay.first - P) / fr);
  const s1 = Math.ceil((lay.last - P) / fr);
  const t0 = Math.floor(c0 / fc);
  const t1 = Math.ceil(c1 / fc);
  const cache = lay.cache!;
  const whole = level.shares.length === 0;
  for (let ti = Math.floor(s0 / ImageTile); ti * ImageTile < s1; ti++)
    for (let tj = Math.floor(t0 / ImageTile); tj * ImageTile < t1; tj++) {
      const key = `${lay.key}|${level.rowLevel}.${level.columnLevel}|${ti}.${tj}|${lay.pyramidOf}`;
      let img = cache.get(key);
      const rows = Math.min(ImageTile, level.rows - ti * ImageTile);
      const cols = Math.min(ImageTile, level.columns - tj * ImageTile);
      if (rows <= 0 || cols <= 0) continue;
      if (!img) {
        img = cache.slot(cols, rows);
        const d = img.data;
        for (let r = 0; r < rows; r++)
          for (let c = 0; c < cols; c++) {
            const o = (ti * ImageTile + r) * level.columns + tj * ImageTile + c;
            const cl = level.classes[o];
            const q = 4 * (r * cols + c);
            if (cl === NoClass) {
              d[q] = d[q + 1] = d[q + 2] = d[q + 3] = 0;
              continue;
            }
            const fill = lay.classColors[cl];
            // opaque by the majority's share: a mixed sample reads lighter (the base level holds no shares: whole)
            const share = whole ? 255 : level.shares[o];
            const a = ((fill & 0xff) * (0.35 + (0.65 * share) / 255)) | 0;
            d[q] = fill >>> 24;
            d[q + 1] = (fill >>> 16) & 0xff;
            d[q + 2] = (fill >>> 8) & 0xff;
            d[q + 3] = a;
          }
        cache.put(key, img);
      }
      const rowStart = ti * ImageTile * fr;
      const colStart = tj * ImageTile * fc;
      const x0 = axis.x(Math.min(colStart, lay.columns));
      const x1 = axis.x(Math.min(colStart + cols * fc, lay.columns));
      const y0 = lay.rowY(P + rowStart);
      const y1 = lay.rowY(P + Math.min(rowStart + rows * fr, R));
      lay.images.push(img);
      lay.ib.push(Math.min(x0, x1), y0, Math.abs(x1 - x0), y1 - y0);
    }
}
