/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The tracks layout on the section layout interface (_spec/README.md "Render", Tracks layout): one strip a section
// (a block in view, once blocks exist), its bands the section's tracks stacked with `trackGap` between them, its frame
// the page (a translation by nothing); its chrome the label column (each track's label, the "n filtered" chip, the
// direction mark). Its viewport is a window of display units across the data box.

import type { Viewport } from '../axis.js';
import type { RenderConfig } from '../config.js';
import { Anchor, type PassName, type PassWriter } from '../draw-lists.js';
import type { Emphasis, Hit } from '../layer-type.js';
import type { SeqTokens } from '../tokens.js';
import type { SectionLayout, TrackLayout } from '../tracks.js';
import { StripPlacement } from './strips.js';
import type { NavContext, PlaceContext, SectionLayoutType, SectionNav, StepDir } from './types.js';

// the anchor a number read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;

/** At least `width / maxScale` wide (within the total), centered, within the axis; a cyclic window may pass the end. */
export function clampWindow(c: NavContext, v: Viewport): Viewport {
  const total = c.total;
  const w = Math.min(Math.max(v.end - v.start, Math.min(c.width / c.maxScale, total)), total);
  if (c.cyclic) {
    const s = ((v.start % total) + total) % total;
    return s === v.start && w === v.end - v.start ? v : { ...v, start: s, end: s + w };
  }
  if (w === v.end - v.start && v.start >= 0 && v.end <= total) return v;
  const center = (v.start + v.end) / 2;
  const s = Math.min(Math.max(center - w / 2, 0), total - w);
  return { ...v, start: s, end: s + w };
}

/** The tracks layout's arithmetic: a window of display units across `width` px. */
export const TracksNav: SectionNav = {
  pan(c, v, dx, dy) {
    const d = (dx * (v.end - v.start)) / c.width;
    const top = Math.max(0, (v.top ?? 0) + dy);
    return clampWindow(c, { ...v, start: v.start - d, end: v.end - d, ...(dy ? { top } : {}) });
  },
  zoom(c, v, factor, at) {
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
    const s = (start + end) / 2 - w / 2;
    return clampWindow(c, { ...v, start: s, end: s + w });
  },
  step(c, offset, dir: StepDir) {
    const last = c.total - 1;
    if (dir === 'home') return 0;
    if (dir === 'end') return last;
    const d = dir === 'left' || dir === 'up' ? -1 : 1;
    const next = offset + d;
    if (c.cyclic) return ((next % c.total) + c.total) % c.total;
    return next < 0 || next > last ? -1 : next;
  },
  clamp: clampWindow,
};

/** Two windows place alike: the same start and end, scroll and bands. */
export function sameWindow(a: Viewport | undefined, b: Viewport | undefined): boolean {
  if (a === b) return true;
  if (!a || !b) return false;
  return a.start === b.start && a.end === b.end && (a.top ?? 0) === (b.top ?? 0) && a.bands === b.bands;
}

export const TracksLayout: SectionLayoutType<StripPlacement> = {
  name: 'tracks',
  fillBands: true,
  place(_section, ctx: PlaceContext, _viewport, prev) {
    const p = prev ?? new StripPlacement();
    const strips = p.strips;
    const { config } = ctx;
    strips.reset(ctx.bands.length);
    const top = ctx.box.y;
    const i = strips.add(ctx.axis, 0, 0, 0, ctx.axis.u0);
    let y = top;
    let bottom = top;
    for (let j = 0; j < ctx.bands.length; j++) {
      if (j) y += config.trackGap;
      strips.bandY[i * strips.bands + j] = y;
      strips.bandH[i * strips.bands + j] = ctx.heights[j];
      y += ctx.heights[j];
      bottom = y;
    }
    strips.top[i] = top;
    strips.height[i] = bottom - top;
    p.contentHeight = bottom - top;
    p.box = { x: config.inset, y: top, width: ctx.width - 2 * config.inset, height: bottom - top };
    ctx.layStrip(p, i);
    return p;
  },
  sameKey: sameWindow,
  nav: TracksNav,
  chrome(s, pass, out, _emphasis, _strip, band) {
    if (pass !== 'labels' || band < 0) return;
    const t = s.tracks[band];
    if (!t || t.box.y >= s.pageHeight) return;
    drawLabel(out, t, s.tokens, s.config);
    if (band === 0 && s.mark) drawMark(out, s, t, s.tokens, s.config);
  },
  hitChrome(s, x, y) {
    for (const t of s.tracks) {
      // the "n filtered" chip in the label column: a click clears the filter
      const f = t.filter;
      if (f && x >= f.x - 2 && x < f.x + f.width + 2 && y >= t.box.y && y < t.box.y + t.box.height)
        return { layer: f.layer, section: s.id, system: '', target: 'filter', count: f.count } satisfies Hit;
    }
    return undefined;
  },
};

/** A track's label in the label column, faded past it. */
function drawLabel(out: PassWriter, t: TrackLayout, tokens: SeqTokens, config: RenderConfig) {
  const g = out.args;
  g[0] = config.inset;
  g[1] = t.labelY;
  g[2] = t.labelWidth;
  const f = t.filter;
  // the label fades before the chip
  g[3] = f ? Math.max(0, f.x - config.inset - 6) : config.labelWidth;
  g[4] = tokens.colors['msui-text-muted'];
  out.textArgs(t.track.label, AnchorStart, config.labelFont);
  if (!f) return;
  g[0] = f.x;
  g[1] = t.labelY;
  g[2] = f.width;
  g[3] = f.width;
  g[4] = tokens.colors['msui-accent'];
  out.textArgs(f.text, AnchorStart, config.labelFont);
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
