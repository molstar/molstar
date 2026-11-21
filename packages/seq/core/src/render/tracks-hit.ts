/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import type { Box, Hit } from './layer-type.js';
import { breakHit } from './marks.js';
import { BlockPlacement } from './layout/blocks.js';
import { PlacePoint, WarpKind } from './layout/strips.js';
import { hitRail } from './rail.js';
import type { SectionLayout, ViewLayout } from './tracks-types.js';
import { BreakSlop, Notch } from './tracks-draw.js';

// The view's hit-test (_spec/README.md "Render", Lanes and Hits): a break's notch, the rails, what a layer draws
// outside its band, the section layout's own marks, then the point placed and the topmost layer's hit in its band.

/** The point `hitView` places (reused: a hover makes no object). */
const Point = new PlacePoint();

/**
 * What is under a point (_spec/README.md "Render", Lanes and Hits): a break within `BreakSlop` px of its notch (above
 * a strip's band column, or a hidden items' tick a layer hits as the break), the section layout's own marks, else the
 * point placed (`Placement.at`: its strip, band and the point in the strip's frame) and the topmost layer's hit in that
 * band. On data the cells and items beside a break win: the hairline has no width of its own.
 */
/** The point a layer's hit reads (`LayerType.hit`): one array, so no number is boxed on the call. */
const HitAt = new Float64Array(2);
/** The page point a placement's `at` reads (`Placement.at`: from an array, never passed). */
const PageAt = new Float64Array(2);

/** `HitAt` set to the placed point (in its strip's frame). */
function pointAt(): Float64Array {
  HitAt[0] = Point.x;
  HitAt[1] = Point.y;
  return HitAt;
}

/** A run's block and count (`runAt`), reused. */
const RunAt = new Int32Array(2);

export function hitView(layout: ViewLayout, x: number, y: number, touch = false, hover?: Hit): Hit | undefined {
  const h = hitLayout(layout, x, y, touch, hover);
  // a run of small blocks (_spec/README.md "Render", Blocks): over its band, what names no item (a ruler's place, its
  // band's ground) is the run's: the block under the point and the run's count, its hover "one of 2,841"; a block too
  // narrow to be read says so
  // the ruler rail's place (outside its window: at the fit its whole body) over a run keeps the rail's target (its
  // gestures) and carries the run's count and state (its card "one of 2,841")
  if (h?.target === 'window' && h.part === 'outside' && h.block !== undefined) {
    const r = runHit(layout, x, y);
    return r ? { ...h, block: r.block, system: r.system, count: r.count, ...(r.state ? { state: r.state } : {}) } : h;
  }
  if (h?.element || (h?.target !== undefined && h.target !== 'run')) return h;
  const r = runHit(layout, x, y);
  if (!r) return h;
  // a read block's position stays the hit's (its data drawn at its place)
  return h?.position !== undefined && h.block === r.block ? { ...r, position: h.position } : r;
}

/** A run's hit at a point: the block under it and the run's count (`state` when the block is not read at this zoom). */
function runHit(layout: ViewLayout, x: number, y: number): Hit | undefined {
  for (const s of layout.sections) {
    const p = s.placement;
    if (!(p instanceof BlockPlacement) || !p.runN || y < p.tracksTop || y >= p.tracksBottom || !p.runAt(x, RunAt))
      continue;
    const block = s.section.blocks?.[RunAt[0]];
    if (!block) continue;
    return {
      layer: s.section.address,
      section: s.id,
      system: block.axis.system,
      block: RunAt[0],
      target: 'run',
      count: RunAt[1],
      ...(p.strips.atX(x) < 0 ? { state: 'not read at this zoom' } : {}),
    };
  }
  return undefined;
}

function hitLayout(layout: ViewLayout, x: number, y: number, touch: boolean, hover: Hit | undefined): Hit | undefined {
  for (const s of layout.sections) {
    const strips = s.placement.strips;
    for (let i = 0; i < strips.n; i++) {
      const col = s.columns[i];
      const breaks = s.stripBreaks[i];
      // a polar strip's notches are its layout's to hit (the seam: `hitChrome`, through the warp's inverse)
      if (!col || !breaks.length || strips.warp[i] !== WarpKind.translate) continue;
      // the page point passed as it came (no arithmetic on it here: a number computed is boxed on each call)
      const notch = breakAt(s, i, col, breaks, x, y);
      if (notch) return notch;
    }
    // the ruler rail's grab zone, over the ruler (an idle rail is no target: a press falls through)
    if (s.rails)
      for (let i = 0; i < s.rails.length; i++) {
        const h = hitRail(s.rails[i], x, y, s.config, touch);
        if (!h) continue;
        // the section's other rails show the same place marked
        for (let k = 0; k < s.rails.length; k++) s.rails[k].markAt = s.rails[i].markAt;
        return h;
      }
    if (s.rowsRails)
      for (let i = 0; i < s.rowsRails.length; i++) {
        const h = hitRail(s.rowsRails[i], x, y, s.config, touch);
        if (h) return h;
      }
    // what a layer draws outside its band (a rows layer's column titles, above the rows), on its first strip
    for (let j = 0; j < s.tracks.length; j++) {
      const lays = s.tracks[j].layers;
      for (let i = 0; i < lays.length; i++) {
        const t = lays[i].type;
        if (!t.hitOutside) continue;
        // in the first strip's frame (a text figure's first line)
        HitAt[0] = strips.n ? x - strips.dx[0] : x;
        HitAt[1] = strips.n ? y - strips.dy[0] : y;
        const h = t.hitOutside(lays[i].lay, HitAt);
        if (h) return h;
      }
    }
    const own = s.layout.hitChrome?.(s, x, y, hover, touch);
    if (own) return own;
    // a `section` layer's marks, on the page, over the strips'
    for (let i = s.sectionLays.length - 1; i >= 0; i--) {
      const l = s.sectionLays[i];
      HitAt[0] = x;
      HitAt[1] = y;
      const h = l.type.hit(l.lay, HitAt);
      if (h?.element) return h;
    }
    PageAt[0] = x;
    PageAt[1] = y;
    if (!s.placement.at(PageAt, Point)) continue;
    // a band's header in the label column (a rows layer's tree and labels): its layers of two axes hit it
    if (Point.kind === 'gutter' && Point.band >= 0) {
      const lays = s.lays[Point.strip * strips.bands + Point.band];
      for (let i = lays ? lays.length - 1 : -1; i >= 0; i--) {
        if (lays![i].type.axes !== 2) continue;
        const h = lays![i].type.hit(lays![i].lay, pointAt());
        if (h) return withBlock(s, Point.strip, h);
      }
    }
    if (Point.kind !== 'data') continue;
    // a layer drawn across other bands of its strip (a frame's outline, `across`) is hit there, over their layers
    for (let j = 0; j < strips.bands; j++) {
      const own = s.lays[Point.strip * strips.bands + j];
      for (let i = own ? own.length - 1 : -1; i >= 0; i--) {
        if (!Array.isArray(own![i].layer.params.across)) continue;
        const h = own![i].type.hit(own![i].lay, pointAt());
        if (h?.element) return h;
      }
    }
    const lays = s.lays[Point.strip * strips.bands + Point.band];
    if (!lays) return undefined;
    // the topmost layer with an element under the point, else the topmost with a position
    let position: Hit | undefined;
    for (let i = lays.length - 1; i >= 0; i--) {
      const l = lays[i];
      const h = l.type.hit(l.lay, pointAt());
      if (h?.element || h?.target === 'break') return withBlock(s, Point.strip, h);
      position ??= h;
    }
    return position && withBlock(s, Point.strip, position);
  }
  return undefined;
}

/**
 * A hit with its block (its index in the section) when the section has blocks (_spec/README.md "Render", Lanes and
 * Hits): written into the layer's hit, no copy (a layer's hit is its lay's, and a lay is one strip's: one block's).
 */
function withBlock(s: SectionLayout, strip: number, h: Hit): Hit {
  if (s.section.blocks) {
    const st = s.placement.strips;
    (h as { block?: number }).block = st.block[strip];
    // a grid over blocks: the strip's y block too
    if (st.ys[strip]) (h as { yBlock?: number }).yBlock = st.yBlock[strip];
  }
  return h;
}

/** The break whose notch is within `BreakSlop` px of a point above the band column (in the strip's frame); else undefined. */
/** A notch's hit at a page point over strip `i`. */
function breakAt(
  s: SectionLayout,
  i: number,
  band: Box,
  b: readonly number[],
  px: number,
  py: number,
): Hit | undefined {
  const strips = s.placement.strips;
  const axis = strips.axes[i];
  const x = px - strips.dx[i];
  const y = py - strips.dy[i];
  if (y < band.y - Notch.height - BreakSlop || y >= band.y) return undefined;
  for (let k = 0; k + 2 < b.length; k += 3)
    if (Math.abs(x - b[k]) <= BreakSlop) return breakHit(s, axis, b[k], b[k + 1], b[k + 2]);
  return undefined;
}
