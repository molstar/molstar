/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tracks layout on the section layout interface (_spec/README.md "Render", Tracks layout): one strip a section
// (a block in view, once blocks exist), its bands the section's tracks stacked with `trackGap` between them, its frame
// the page (a translation by nothing); its chrome the label column (each track's label and header line, the section's
// header line, the direction mark). Its viewport is a window of display units across the data box.

import type { Viewport } from '../axis.js';
import type { RenderConfig } from '../config.js';
import { Anchor, type PassName, type PassWriter } from '../draw-lists.js';
import type { Emphasis } from '../layer-type.js';
import type { SeqTokens } from '../tokens.js';
import type { SectionLayout, TrackLayout } from '../tracks.js';
import { BlockPlacement, blockChrome, placeBlocks, rowPan, rowZoom } from './blocks.js';
import { capScroll, headOf, StripPlacement } from './strips.js';
import { drawSectionHeader, drawTrackHeader, hitHeaders } from '../headers-draw.js';
import { drawGroupHeaders, hitGroupHeader } from '../groups-draw.js';
import type { NavContext, PlaceContext, SectionLayoutType, SectionNav, StepDir } from './types.js';

// the anchor a number read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;

/** At least `width / maxScale` wide (within the total), centered, within the axis; a cyclic window may pass the end. */
export function clampWindow(c: NavContext, v: Viewport): Viewport {
  const t = clampTop(c, v.top);
  if (t !== v.top) v = { ...v, top: t };
  const total = c.total;
  const w = Math.min(Math.max(v.end - v.start, Math.min(c.width / c.maxScale, total)), total);
  if (c.cyclic) {
    // as wide as the axis: the fit (a whole turn from elsewhere would show the origin mid-window, a follower over it)
    if (w >= total) return v.start === 0 && v.end === total ? v : { ...v, start: 0, end: total };
    const s = ((v.start % total) + total) % total;
    return s === v.start && w === v.end - v.start ? v : { ...v, start: s, end: s + w };
  }
  if (w === v.end - v.start && v.start >= 0 && v.end <= total) return v;
  const center = (v.start + v.end) / 2;
  const s = Math.min(Math.max(center - w / 2, 0), total - w);
  return { ...v, start: s, end: s + w };
}

/**
 * A section's scroll within its content (_spec/README.md "Render", Section scroll): `[0, content − box]` where it
 * scrolls, 0 where it shows all of it; as given with no placement.
 */
function clampTop(c: NavContext, top: number | undefined): number | undefined {
  const p = c.placement;
  if (top === undefined || !p) return top;
  const max = p.clip ? Math.max(0, p.contentHeight - p.box.height) : 0;
  const t = Math.max(0, Math.min(max, top));
  return t === 0 && top === 0 ? top : t;
}

/** The tracks layout's arithmetic: a window of display units across `width` px. */
export const TracksNav: SectionNav = {
  pan(c, v, dx, dy, x) {
    const top = Math.max(0, (v.top ?? 0) + dy);
    // a row of blocks pans through its inverse: what is under the pointer moves with it over chips and gaps
    const row = rowOf(c);
    const moved = row ? rowPan(row, v, dx, x) : undefined;
    if (moved) return clampWindow(c, { ...moved, ...(dy ? { top } : {}) });
    const d = (dx * (v.end - v.start)) / c.width;
    return clampWindow(c, { ...v, start: v.start + d, end: v.end + d, ...(dy ? { top } : {}) });
  },
  zoom(c, v, factor, at) {
    // a row of blocks keeps the anchor's px (chips entering or leaving the view move it not)
    const row = rowOf(c);
    if (row) return clampWindow(c, rowZoom(row, c, v, factor, at));
    const start = at - (at - v.start) / factor;
    const end = at + (v.end - at) / factor;
    return clampWindow(c, { ...v, start, end });
  },
  fit(c, v) {
    return { ...v, start: 0, end: c.total };
  },
  reveal(c, v, start, end) {
    if (start >= v.start && end <= v.end) return v;
    const w = v.end - v.start;
    if (end - start > w) return clampWindow(c, { ...v, start, end });
    // the nearer edge: what is left of the window comes in at its left edge, what is right at its right
    const d = start < v.start ? start - v.start : end - v.end;
    return clampWindow(c, { ...v, start: v.start + d, end: v.end + d });
  },
  step(c, offset, dir: StepDir) {
    const last = c.total - 1;
    if (dir === 'home') return 0;
    if (dir === 'end') return last;
    // one line: up and down leave it
    if (dir === 'up' || dir === 'down') return -1;
    const d = dir === 'left' ? -1 : 1;
    const next = offset + d;
    if (c.cyclic) return ((next % c.total) + c.total) % c.total;
    return next < 0 || next > last ? -1 : next;
  },
  clamp: clampWindow,
};

/** A row of blocks' arithmetic, when the context's placement is one (a stack is linear). */
function rowOf(c: NavContext) {
  const p = c.placement;
  return p instanceof BlockPlacement && !p.stacked && !c.cyclic && p.row.n > 0 ? p.row : undefined;
}

/** Two windows place alike: the same start and end, scroll and bands (by value). */
export function sameWindow(a: Viewport | undefined, b: Viewport | undefined): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  return (
    a.start === b.start &&
    a.end === b.end &&
    (a.top ?? 0) === (b.top ?? 0) &&
    a.groups === b.groups &&
    sameBands(a.bands, b.bands)
  );
}

/** Two viewports' bands alike: each band's scroll and row zoom. */
export function sameBands(a: Viewport['bands'], b: Viewport['bands']): boolean {
  if (a === b) return true;
  const ka = a ? Object.keys(a) : [];
  const kb = b ? Object.keys(b) : [];
  if (ka.length !== kb.length) return false;
  for (const k of ka) {
    const x = a![k];
    const y = b?.[k];
    if (!y || x.bandTop !== y.bandTop || x.rowScale !== y.rowScale || x.folds !== y.folds) return false;
  }
  return true;
}

export const TracksLayout: SectionLayoutType<StripPlacement> = {
  name: 'tracks',
  fillBands: true,
  // a cyclic axis's window passes its origin: drawn as its two runs (the axis's pieces twice, `layoutView`)
  wraps: true,
  place(section, ctx: PlaceContext, viewport, prev) {
    if (section.blocks)
      return placeBlocks(section, ctx, prev instanceof BlockPlacement ? prev : undefined, viewport?.top);
    const p = prev && !(prev instanceof BlockPlacement) ? prev : new StripPlacement();
    const strips = p.strips;
    const { config } = ctx;
    strips.reset(ctx.bands.length);
    const top = ctx.box.y;
    const i = strips.add(ctx.axis, 0, 0, 0, ctx.axis.u0);
    let y = top;
    let bottom = top;
    const pads = ctx.pads;
    const leads = ctx.leads;
    for (let j = 0; j < ctx.bands.length; j++) {
      if (j) y += config.trackGap;
      if (leads) y += leads[j];
      const k = i * strips.bands + j;
      // a band inside its well: its padding above and below (_spec/README.md "Render", Tracks layout)
      const pt = pads ? pads[2 * j] : 0;
      const pb = pads ? pads[2 * j + 1] : 0;
      strips.padT[k] = pt;
      strips.padB[k] = pb;
      strips.bandY[k] = y + pt;
      strips.bandH[k] = ctx.heights[j];
      y += pt + ctx.heights[j] + pb;
      bottom = y;
    }
    // past its box its tracks below the ruler scroll (_spec/README.md "Render", Section scroll); an export shows all
    const content = bottom - top;
    const fixed = headOf(ctx.bands, ctx.heights);
    const head = fixed > 0 ? strips.bandY[fixed - 1] + strips.bandH[fixed - 1] + strips.padB[fixed - 1] : top;
    const boxH = capScroll(p, ctx, section.params.height, top, content, fixed, head, false, viewport?.top);
    strips.top[i] = top;
    strips.height[i] = p.clip ? boxH : content;
    p.contentHeight = content;
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: boxH };
    ctx.layStrip(p, i);
    return p;
  },
  sameKey: sameWindow,
  nav: TracksNav,
  chrome(s, pass, out, emphasis, strip, band) {
    const p = s.placement;
    if (p instanceof BlockPlacement) {
      blockChrome(s, p, pass, out, strip, band);
      // the label column once a row: a stack's chips share one
      if (strip > 0 && (!p.stacked || !p.rowStart(strip))) return;
      // a band left out of a stack's row has no label
      if (strip >= 0 && band >= 0 && p.strips.bandH[strip * p.strips.bands + band] <= 0) return;
    }
    // the section's header line (_spec/README.md "Render", Track headers), once before its strips
    if (band < 0) {
      if (strip < 0 && (pass === 'labels' || pass === 'tip')) drawSectionHeader(out, s, pass, emphasis);
      if (strip < 0 && pass === 'labels' && s.groups) drawGroupHeaders(out, s);
      return;
    }
    if (pass !== 'labels' && pass !== 'tip') return;
    const t = s.tracks[band];
    if (!t || t.box.y >= s.pageHeight) return;
    // a stack's later blocks: each its row's labels at its own place
    const dy = strip > 0 && p instanceof BlockPlacement ? p.strips.bandY[strip * p.strips.bands + band] - t.box.y : 0;
    if (pass === 'labels') {
      drawLabel(out, t, s.tokens, s.config, dy);
      if (s.empties[band] === 'note') drawEmptyNote(out, t, s.tokens, s.config, dy, s.note);
      if (band === 0 && s.mark && strip <= 0) drawMark(out, s, t, s.tokens, s.config);
    }
    // its header's states and, when they show, its actions (in a stack's first row)
    drawTrackHeader(out, s, t, pass, emphasis, dy);
  },
  hitChrome(s, x, y, hover, touch) {
    return hitHeaders(s, x, y, hover, touch) ?? (s.groups ? hitGroupHeader(s, x, y) : undefined);
  },
};

/** A track's label in the label column, faded past it (`dy` below its first row: a stack's later blocks). */
function drawLabel(out: PassWriter, t: TrackLayout, tokens: SeqTokens, config: RenderConfig, dy = 0) {
  const g = out.args;
  // a group's tracks' labels indented (_spec/README.md "Render", Track groups)
  const indent = t.indent ?? 0;
  g[0] = config.inset + indent;
  g[1] = t.labelY + dy;
  g[2] = t.labelWidth;
  // the label fades before the states beside it (states on a line of their own leave it whole)
  g[3] = (t.header ? t.header.label.width : config.labelWidth) - indent;
  g[4] = tokens.colors['msui-text-muted'];
  // an empty label (`label: ''`) draws no text
  if (t.labelText) out.textArgs(t.labelText, AnchorStart, config.labelFont);
}

/** What a track with no items says in its row, quiet: the host's "none found" stays visible (_spec/README.md "Render", Tracks layout). */
export const EmptyNote = 'no items';

/** A track's note, at the start of its data box, centred on its row. */
function drawEmptyNote(
  out: PassWriter,
  t: TrackLayout,
  tokens: SeqTokens,
  config: RenderConfig,
  dy = 0,
  note = EmptyNote,
) {
  const g = out.args;
  g[0] = t.box.x + 8;
  g[1] = t.labelY + dy;
  g[2] = Math.max(0, t.box.width - 16);
  g[3] = Math.max(0, t.box.width - 16);
  g[4] = tokens.colors['msui-text-muted'];
  out.textArgs(note, AnchorStart, config.itemFont);
}

/** The direction mark after the first track's label, muted, faded at the label column's edge. */
function drawMark(out: PassWriter, s: SectionLayout, t: TrackLayout, tokens: SeqTokens, config: RenderConfig) {
  const m = s.mark!;
  const g = out.args;
  g[0] = m.x;
  g[1] = t.labelY;
  g[2] = m.width;
  g[3] = Math.max(0, config.inset + config.labelWidth - m.x);
  g[4] = tokens.colors['msui-text-muted'];
  out.textArgs(m.text, AnchorStart, config.labelFont);
}

export type { Emphasis, PassName };
