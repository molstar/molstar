/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Developed in collaboration with Anyfold <https://anyfold.ai/>
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The composer's entries (draw/composer.ts): an entry's parts in the clip (its header, data columns, data and
// gutter), a track's layers through their registry entries (`draw`, `tiny` below its minimum height, else bars; a
// throw omits the layer, reported once), a split row's guide, the rows' band in one `rows` call, and a moved lazy
// band's cut rows. Per frame nothing is allocated (RN1's draw pattern).

import { Lanes } from '@molstar/seq-core/format';
import { AxisGeometry } from '../axis.js';
import type { ViewModel } from '../graph.js';
import { EntryKind, type SectionLayout, SectionLayouts, SPLIT_GUIDE_PX } from '../layout.js';
import type { ViewportState } from '../viewport.js';
import { type RenderBackend, SkinToken } from './backend.js';
import {
  banded,
  bands,
  boxOf,
  breaks,
  type ComposerFrame,
  fail,
  gutterWidth,
  onRow,
  open,
  Part,
  place,
  prime,
  release,
  setClip,
  stateOf,
} from './composer-frame.js';
import { columns, gutter, header, rowsHeader, trackLabel } from './composer-headers.js';
import { drawGroupLabels } from './numbers.js';

/** A split row's guide: the grid tone at this alpha (0..255), well under a row border's full strength. */
const SPLIT_GUIDE_ALPHA = 150;

// Enum-like members read once (a read per frame allocated under the test runner's module wrapper).
const WRAPPED = AxisGeometry.wrapped;
const LABELS = EntryKind.labels;
const TRACK = EntryKind.track;
const TRACK_LABEL = EntryKind.trackLabel;
const LETTERS = 'letters';
const PART_DATA = Part.data;
const PART_GUTTER = Part.gutter;
const PART_HEADER = Part.header;
const TOKEN_GRID = SkinToken.grid;
const TOKEN_INK_FAINT = SkinToken.inkFaint;

/** The parts of entry e inside the clip. */
export function entry(f: ComposerFrame, l: SectionLayout, e: number, parts: number) {
  const v = f.v;
  const b = f.box;
  if (!boxOf(f, l, e)) return;
  f.stats.entries++;
  if (banded(f, l, e)) return band(f, l, parts);
  if (parts & PART_HEADER) header(f, l, e);
  if (parts & PART_HEADER && l.rows.layer >= 0 && l.entryTrack[e] === l.rows.track && e >= l.rows.entry) {
    f.rg[0] = b[0] - l.entryRow[e] * b[1];
    f.rg[1] = b[1];
    columns(f, l, l.entryRow[e], l.entryRow[e]);
  }
  const kind = l.entryKind[e];
  if (kind === LABELS && parts & PART_DATA && open(f, l.dataLeft, l.dataRight, parts)) {
    v.top = b[0];
    v.height = b[1];
    place(f, l, e);
    drawGroupLabels(v);
    v.backend.popClip();
  }
  if (kind === TRACK_LABEL && parts & PART_DATA && open(f, l.dataLeft, l.dataRight, parts)) {
    trackLabel(f, l, e);
    v.backend.popClip();
  }
  if (kind !== TRACK) return;
  if (parts & PART_DATA) data(f, l, e, parts);
  if (parts & PART_GUTTER) gutter(f, l, e);
}
/** The track's layers (an expanded split track's row: its layer; its label's entry: the others), then its border. */
function data(f: ComposerFrame, l: SectionLayout, e: number, parts: number) {
  // A chip's entry is its chip alone (seqspec:viewer#the-text-view): its clip, and a clear, never reach the chips
  // beside it on its row.
  let x0 = l.dataLeft;
  let x1 = l.dataRight;
  if (l.geometry === WRAPPED && l.entryRow[e] >= 0) {
    const line = l.entryRow[e];
    const run = f.v.runs[l.lines.block[line]];
    if (run?.chip) {
      x0 += l.lines.left[line];
      x1 = x0 + run.advance;
    }
  }
  if (!open(f, x0, x1, parts)) return;
  const v = f.v;
  const b = f.box;
  const rv = (f.model as ViewModel).resolved;
  v.top = b[0];
  v.height = b[4];
  place(f, l, e);
  // Wrapped, with the line's first letters: what the line's other tracks draw on its text, under the letters (over
  // them where the letters sit on cells, `wrappedCells`: a highlight is translucent).
  const text = l.geometry === WRAPPED && lettersEntry(f, l, e);
  const under = !v.props.wrappedCells;
  if (text && under) behindText(f, l, e);
  // A split row's guide, under its data.
  if (b[4] === b[1] && b[1] * v.dpr >= SPLIT_GUIDE_PX && SectionLayouts.splitGuide(l, e)) guide(f, b);
  const t = l.entryTrack[e];
  const k0 = l.entryLayer[e];
  const one = SectionLayouts.splitRow(l, e) >= 0;
  if (l.geometry === WRAPPED && t === l.textRows.track) textRows(f, l, e);
  else
    for (let k = k0; k < rv.layers.count && rv.layers.track[k] === t && !(one && k > k0); k++)
      if (l.entryRow[e] !== -2 || k !== rv.tracks.split[t]) layer(f, k, l, e);
  if (text && !under) behindText(f, l, e);
  // The row's border, inside its box below its content (thin split rows, a grid, share one: the track's bottom).
  if (b[4] < b[1]) {
    const r = f.rects;
    const w = Math.max(1, Math.round(v.dpr)) | 0;
    r.reset();
    r.begin(TOKEN_GRID);
    const o = f.opened;
    r.push(o[0], (Math.round((b[0] + b[1]) * v.dpr) - w) | 0, o[1] - o[0], w);
    v.backend.rects(r);
  }
  v.backend.popClip();
}
/** A split row's guide: a hairline at the row's bottom edge in the grid tone, over the opened clip's width. */
function guide(f: ComposerFrame, b: Float64Array) {
  const r = f.rects;
  const w = Math.max(1, Math.round(f.v.dpr)) | 0;
  r.reset();
  r.begin(TOKEN_GRID, SPLIT_GUIDE_ALPHA);
  const o = f.opened;
  r.push(o[0], (Math.round((b[0] + b[1]) * f.v.dpr) - w) | 0, o[1] - o[0], w);
  f.v.backend.rects(r);
}
/** One layer through its entry: `draw`, `tiny` below its minimum height, else bars; a throw omits it. */
function layer(f: ComposerFrame, k: number, l: SectionLayout, e: number) {
  const r = f.renderers[k];
  if (!r || f.failed.has(k)) return;
  const v = f.v;
  const s = stateOf(f, k, l, e);
  const p = (f.model as ViewModel).resolved.layers.params[k];
  prime(f, k, l, e);
  try {
    if (!r.draw) bare(f, k, l);
    else if (v.ring && !r.geometries.includes('circular')) ringless(f, k);
    else if (l.geometry === WRAPPED || v.height >= (r.minHeight ?? 0)) {
      const more = v.more;
      v.more = false;
      r.draw(v, s, p);
      if (v.more) f.pend[k] = 1;
      v.more = v.more || more;
    } else if (r.tiny) r.tiny(v, s, p);
    else bars(f);
  } catch (err) {
    fail(f, k, err);
  } finally {
    release(f);
  }
}
/**
 * Wrapped, the rows drawn as text (seqspec:viewer#rows-in-text): the rows layer over the line's band, then each other
 * layer of the track on the text row of the row it targets (`v.row`, its lane offset), its `behindText` alone (a range's
 * `box` or `highlight`); the rows in view only. A throw omits the layer.
 */
function textRows(f: ComposerFrame, l: SectionLayout, e: number) {
  const v = f.v;
  const tr = l.textRows;
  const L = (f.model as ViewModel).resolved.layers;
  const lanes = (f.model as ViewModel).resolved.lanes;
  layer(f, tr.layer, l, e);
  const top = v.top;
  const height = v.height;
  const h = tr.height;
  const d0 = Math.max(0, Math.floor((f.box[2] - top) / h)) | 0;
  const d1 = Math.min(tr.count, Math.ceil((f.box[3] - top) / h)) | 0;
  const t = l.entryTrack[e];
  for (let k = Math.max(0, l.entryLayer[e]); k < L.count && L.track[k] === t; k++) {
    const r = f.renderers[k];
    if (k === tr.layer || !r?.behindText || f.failed.has(k)) continue;
    prime(f, k, l, e);
    try {
      for (let d = d0; d < d1; d++) {
        v.top = top + d * h;
        v.height = h;
        v.row = d;
        v.dataRow = Lanes.rowOf(lanes, tr.laneStart + d);
        r.behindText(v, stateOf(f, k, l, e), L.params[k]);
      }
    } catch (err) {
      fail(f, k, err);
    } finally {
      release(f);
    }
  }
  v.top = top;
  v.height = height;
  v.row = v.dataRow = -1;
}
/** Whether entry e's track holds a `letters` layer: the text's band on a wrapped line. */
function lettersEntry(f: ComposerFrame, l: SectionLayout, e: number): boolean {
  const L = (f.model as ViewModel).resolved.layers;
  const t = l.entryTrack[e];
  for (let k = Math.max(0, l.entryLayer[e]); k < L.count && L.track[k] === t; k++)
    if (f.renderers[k]?.name === LETTERS) return true;
  return false;
}
/**
 * Wrapped, letters entry e: when it is its line's first, each layer of the line's tracks with `behindText` draws on the
 * text's band (`v.top`, `v.height` as placed), with its own state and themes; a throw omits the layer.
 */
function behindText(f: ComposerFrame, l: SectionLayout, e: number) {
  const line = l.entryRow[e];
  const top = l.lines.top[line];
  const bottom = l.lines.bottom[line];
  for (let q = e - 1; q >= l.pinnedCount && l.entryTop[q] >= top; q--)
    if (l.entryKind[q] === TRACK && l.entryRow[q] === line && lettersEntry(f, l, q)) return;
  const v = f.v;
  const L = (f.model as ViewModel).resolved.layers;
  for (let q = e + 1; q < l.entryCount && l.entryTop[q] < bottom; q++) {
    if (l.entryKind[q] !== TRACK || l.entryRow[q] !== line) continue;
    const t = l.entryTrack[q];
    for (let k = l.entryLayer[q]; k < L.count && L.track[k] === t; k++) {
      const r = f.renderers[k];
      if (!r?.behindText || f.failed.has(k)) continue;
      prime(f, k, l, q);
      try {
        r.behindText(v, stateOf(f, k, l, q), L.params[k]);
      } catch (err) {
        fail(f, k, err);
      } finally {
        release(f);
      }
    }
  }
}
/** A renderer that declares the geometry and has no `draw`: reported once, drawn as bars. */
function bare(f: ComposerFrame, k: number, l: SectionLayout) {
  if (!f.undrawn.has(k)) {
    f.undrawn.add(k);
    const geometry = l.geometry === WRAPPED ? 'wrapped' : 'linear';
    f.ctx.report('svs.render.fallback', `layer "${f.refs[k]}" declares ${geometry} and has no draw: drawn as bars`);
  }
  bars(f);
}
/**
 * On rings, a renderer that does not declare `circular` (it draws images, which a ring does not show): its extents as
 * bars, and the section says so once per layer (`svs.render.fallback`).
 */
function ringless(f: ComposerFrame, k: number) {
  if (!f.unringed.has(k)) {
    f.unringed.add(k);
    const m = f.model as ViewModel;
    const type = m.set.strings.get(m.resolved.layers.type[k]);
    f.ctx.report(
      'svs.render.fallback',
      `layer "${f.refs[k]}" (${type}) has no ring form in section "${f.section}": drawn as bars`,
    );
  }
  bars(f);
}
/** Below a renderer's minimum height, without `tiny`: its blocks' shown extents as bars. */
function bars(f: ComposerFrame) {
  const v = f.v;
  const r = f.rects;
  const y = Math.round((v.top + v.height * 0.25) * v.dpr) | 0;
  const h = Math.max(1, Math.round(v.height * 0.5 * v.dpr)) | 0;
  r.reset();
  r.begin(TOKEN_INK_FAINT);
  for (let b = 0; b < v.visLo.length; b++) {
    if ((v.onlyBlock >= 0 && b !== v.onlyBlock) || v.visHi[b] <= v.visLo[b]) continue;
    const x = v.xd(b, v.visLo[b]);
    r.push(x, y, v.xd(b, v.visHi[b]) - x, h);
  }
  if (r.n > 0) v.backend.rects(r);
}
/**
 * The band may be copied: it draws no tree (paths, which a strip's clip rasterizes differently). Its tiles are drawn
 * at their whole scale and a pending one's coarse level in its place alone, so a strip draws what a full draw does,
 * regions in view pending or not; an arrival or a labels chunk redraws what changed.
 */
export function bandCopies(f: ComposerFrame, l: SectionLayout): boolean {
  const hd = f.renderers[l.rows.layer]?.header;
  const m = f.model as ViewModel;
  return !hd?.draw || hd.inset(m.prepared[l.rows.layer].state, m.resolved.layers.params[l.rows.layer]) <= 0;
}
/**
 * A moved lazy band's rows a scrolled region's edge cuts, each whole under the clip a full draw has (not the band);
 * a dense band draws no letters, and its images are cut alike.
 */
export function cutRows(
  f: ComposerFrame,
  be: RenderBackend,
  l: SectionLayout,
  s: ViewportState,
  W: number,
  y0: number,
  y1: number,
) {
  const r = l.rows;
  const dpr = f.v.dpr;
  if (r.height < f.v.props.denseRowHeight) return;
  const top = l.scrolledTop + l.entryTop[r.entry] - s.scrollY;
  for (let side = 0; side < 2; side++) {
    const edge = side === 0 ? l.scrolledTop : l.height;
    const d = Math.floor((edge - top) / r.height);
    if (d < 0 || d >= r.count || top + d * r.height >= edge) continue;
    const a = Math.max(y0, Math.floor((top + d * r.height) * dpr));
    const b = Math.min(y1, Math.ceil((top + (d + 1) * r.height) * dpr));
    if (b <= a) continue;
    be.clear(0, a, W, b - a);
    setClip(f, 0, W, a, b);
    bands(f, l);
    f.dense = false;
    entry(f, l, r.entry, PART_HEADER | PART_DATA);
    breaks(f, l);
  }
}
/** The rows' band in the clip: the display rows in it in one `rows` call (once per region), and their header. */
function band(f: ComposerFrame, l: SectionLayout, parts: number) {
  const r = l.rows;
  if (f.dense) return;
  f.dense = true;
  const v = f.v;
  const h = r.height;
  const top0 = l.scrolledTop + l.entryTop[r.entry] - v.scrollY;
  const y0 = Math.max(l.scrolledTop, f.clip[2] / v.dpr);
  const y1 = Math.min(l.height, f.clip[3] / v.dpr);
  const d0 = Math.max(0, Math.floor((y0 - top0) / h));
  const d1 = Math.min(r.count - 1, Math.ceil((y1 - top0) / h) - 1);
  if (d1 < d0) return;
  const b = f.box;
  b[0] = top0 + d0 * h;
  b[1] = (d1 - d0 + 1) * h;
  b[2] = y0;
  b[3] = y1;
  f.rg[0] = top0;
  f.rg[1] = h;
  // The rows' labels within the band's visible part: a row the scrolled region's top cuts draws no label above it.
  if (parts & PART_HEADER && open(f, 0, Math.min(l.headerWidth, l.dataLeft - gutterWidth(f, l)), 0)) {
    rowsHeader(f, l, d0, d1);
    v.backend.popClip();
  }
  if (parts & PART_HEADER) columns(f, l, d0, d1);
  const k = r.layer;
  const rend = f.renderers[k];
  if (!(parts & PART_DATA) || !rend?.rows || f.failed.has(k) || !open(f, l.dataLeft, l.dataRight, parts)) return;
  const m = f.model as ViewModel;
  prime(f, k, l, r.entry);
  v.top = top0 + d0 * h;
  v.height = h;
  v.row = v.dataRow = v.onlyBlock = -1;
  f.stats.rowsCalls++;
  try {
    rend.rows(v, m.prepared[k].state, m.resolved.layers.params[k], d0, d1);
  } catch (err) {
    fail(f, k, err);
  } finally {
    release(f);
  }
  // The track's other layers on each row's lane, where rows are tall enough to read them (lazy rows).
  const rv = m.resolved;
  for (let q = k + 1; h >= v.props.denseRowHeight && q < rv.layers.count && rv.layers.track[q] === r.track; q++)
    for (let d = d0; d <= d1; d++) {
      v.top = top0 + d * h;
      v.height = h;
      onRow(f, l, d);
      layer(f, q, l, r.entry);
    }
  v.backend.popClip();
}
