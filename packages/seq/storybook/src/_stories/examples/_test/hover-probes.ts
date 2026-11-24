/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import {
  DrawLists,
  type Frame,
  type GlyphList,
  type ImageList,
  type PassList,
  type PathList,
  PathCommandPoints,
  type RectList,
  type SectionLayout,
  walkPass,
  WarpKind,
} from '@molstar/seqx-core/render';

// The hover contract's sweep (_spec/README.md "Render", Emphasis, The hover contract): points on the marks each layer
// instance drew, read from its own data pass (drawn again alone into scratch lists, so a mark is a layer's by
// construction, never guessed from the frame's merged lists), carried to the page through its strip's warp; the
// pointer moved to each, and what the viewer answers checked.

/** A point on a drawn mark of a layer: its page point, the section it is in and the layer's key and type. */
export interface Probe {
  readonly x: number;
  readonly y: number;
  /** Its section's address, or a bridge's id. */
  readonly section: string;
  readonly layer: string;
  readonly type: string;
  /** What the mark is: `rect`, `path`, `letter`, `image`. */
  readonly mark: string;
  /** Its band's top and foot on the page, in a translated strip (absent on a ring, a section's own layer, a bridge). */
  readonly band?: { readonly y0: number; readonly y1: number };
}

type Point = { x: number; y: number; mark: string };
export type Bounds = { x0: number; y0: number; x1: number; y1: number };

const scratch = new DrawLists();

/** Points on a pass's marks, in its frame: a rect's centre, a path's middle point on its curve, a letter's, an image's. */
function marksOf(p: PassList, out: Point[]) {
  walkPass(p, {
    rect(l: RectList, i: number) {
      if (l.w[i] < 0.5 || l.h[i] < 0.5) return;
      out.push({ x: l.x[i] + l.w[i] / 2, y: l.y[i] + l.h[i] / 2, mark: 'rect' });
    },
    path(l: PathList, i: number) {
      // the end point of its middle command: on the curve (a cubic's control points are not)
      const c0 = l.cmdStart[i];
      const n = l.cmdCount[i];
      if (!n) return;
      let pt = l.pointStart[i];
      const mid = c0 + Math.floor(n / 2);
      for (let c = c0; c < mid; c++) pt += PathCommandPoints[l.commands[c]];
      const k = l.commands[mid] === 3 ? pt - 1 : pt + PathCommandPoints[l.commands[mid]] - 1;
      if (k < 0) return;
      // a pixel in from the extent drawn so far (an area's baseline point lies on its band's foot, outside the band)
      let x0 = Infinity;
      let y0 = Infinity;
      let x1 = -Infinity;
      let y1 = -Infinity;
      for (let q = l.pointStart[i]; q <= k; q++) {
        x0 = Math.min(x0, l.points[2 * q]);
        x1 = Math.max(x1, l.points[2 * q]);
        y0 = Math.min(y0, l.points[2 * q + 1]);
        y1 = Math.max(y1, l.points[2 * q + 1]);
      }
      const x = l.points[2 * k];
      const y = l.points[2 * k + 1];
      out.push({
        x: Math.min(Math.max(x, x0 + 1), x1 - 1) || x,
        y: Math.min(Math.max(y, y0 + 1), y1 - 1) || y,
        mark: 'path',
      });
    },
    glyphs(l: GlyphList, i: number) {
      const j = l.start[i] + Math.floor(l.count[i] / 2);
      if (!l.count[i]) return;
      // a letter's centre, a little above its baseline
      const y = l.placed[i] ? l.ly[j] : l.y[i];
      out.push({ x: l.x[j], y: y - 3, mark: 'letter' });
    },
    text() {},
    image(l: ImageList, i: number) {
      // the pixels it paints (a density strip's empty bins are clear): a few spread over them, at their centres
      const im = l.images[i];
      const lit: number[] = [];
      for (let k = 0; k < im.width * im.height; k++) if (im.data[4 * k + 3] > 0) lit.push(k);
      for (const k of spread(lit, 3)) {
        const px = k % im.width;
        const py = Math.floor(k / im.width);
        out.push({
          x: l.x[i] + ((px + 0.5) * l.w[i]) / im.width,
          y: l.y[i] + ((py + 0.5) * l.h[i]) / im.height,
          mark: 'image',
        });
      }
    },
    clip() {},
    unclip() {},
  });
}

/** Up to `per` of a list, spread over it (its first, its last and evenly between). */
function spread<T>(xs: readonly T[], per: number): T[] {
  if (xs.length <= per) return [...xs];
  const out: T[] = [];
  for (let k = 0; k < per; k++) out.push(xs[Math.round((k * (xs.length - 1)) / (per - 1))]);
  return out;
}

/** Whether a page point lies in a section's visible box. */
function inBox(s: SectionLayout, x: number, y: number): boolean {
  const b = s.placement.box;
  return x >= b.x && x < b.x + b.width && y >= b.y && y < b.y + b.height;
}

/**
 * The probes of a frame: for every layer instance (a strip's band's layers, a section's own layers, a bridge's
 * ribbons), up to `per` points on the marks its data pass drew, inside its section's visible box (in its band).
 */
export function probesOf(frame: Frame, per = 4): Probe[] {
  const out: Probe[] = [];
  const pass = scratch.passes.data;
  const w = scratch.writer('data');
  const pts: Point[] = [];
  const take = (
    section: string,
    type: string,
    layer: string,
    map: (x: number, y: number) => [number, number],
    inside: (x: number, y: number) => boolean,
    band?: { y0: number; y1: number },
  ) => {
    pts.length = 0;
    marksOf(pass, pts);
    const placed: Probe[] = [];
    for (const p of pts) {
      const [x, y] = map(p.x, p.y);
      if (!Number.isFinite(x) || !Number.isFinite(y) || !inside(x, y)) continue;
      placed.push({ x, y, section, layer, type, mark: p.mark, ...(band ? { band } : {}) });
    }
    out.push(...spread(placed, per));
  };
  for (const s of frame.layout.sections) {
    const st = s.placement.strips;
    const nb = st.bands;
    const inSection = (x: number, y: number) => inBox(s, x, y);
    for (let i = 0; i < st.n; i++) {
      const polar = st.warp[i] === WarpKind.polar ? s.placement.warpOf?.(i) : undefined;
      const dx = st.dx[i];
      const dy = st.dy[i];
      const map = polar
        ? (x: number, y: number): [number, number] => {
            polar.map(x, y);
            return [polar.px, polar.py];
          }
        : (x: number, y: number): [number, number] => [x + dx, y + dy];
      for (let j = 0; j < nb; j++) {
        // a mark is seen inside its band alone (the data is clipped to its well, a block's to its strip's box; hit
        // tests read the band)
        const by = st.bandY[i * nb + j];
        const bh = st.bandH[i * nb + j];
        const box = s.boxes[i * nb + j];
        const inBand = (x: number, y: number): [number, number] =>
          y >= by && y < by + bh && x >= box.x && x < box.x + box.width ? map(x, y) : [Number.NaN, Number.NaN];
        for (const l of s.lays[i * nb + j] ?? []) {
          if (l.zoomIn) continue;
          scratch.reset(['data']);
          l.type.draw(l.lay, 'data', w);
          const band = polar || l.layer.params.across ? undefined : { y0: by + dy, y1: by + bh + dy };
          take(s.section.address, l.type.name, l.key, l.layer.params.across ? map : inBand, inSection, band);
        }
      }
    }
    for (const l of s.sectionLays) {
      scratch.reset(['data']);
      l.type.draw(l.lay, 'data', w);
      take(s.section.address, l.type.name, l.key, (x, y) => [x, y], inSection);
    }
  }
  const bridges = frame.layout.bridges;
  if (bridges) {
    scratch.reset(['data']);
    bridges.draw('data', w, undefined);
    for (const b of bridges.list)
      take(
        b.id,
        'bridge',
        b.id,
        (x, y) => [x, y],
        (x, y) => bridges.at(x, y) === b.index,
      );
  }
  return out;
}

/** The boxes of a pass's primitives on the page (rects, paths' extents, images): where its emphasis was drawn. */
export function passBoxes(p: PassList): Bounds[] {
  const out: Bounds[] = [];
  walkPass(p, {
    rect(l: RectList, i: number) {
      out.push({ x0: l.x[i], y0: l.y[i], x1: l.x[i] + l.w[i], y1: l.y[i] + l.h[i] });
    },
    path(l: PathList, i: number) {
      let pt = l.pointStart[i];
      let n = 0;
      for (let c = l.cmdStart[i]; c < l.cmdStart[i] + l.cmdCount[i]; c++) n += PathCommandPoints[l.commands[c]];
      let x0 = Infinity;
      let y0 = Infinity;
      let x1 = -Infinity;
      let y1 = -Infinity;
      for (const end = pt + n; pt < end; pt++) {
        x0 = Math.min(x0, l.points[2 * pt]);
        x1 = Math.max(x1, l.points[2 * pt]);
        y0 = Math.min(y0, l.points[2 * pt + 1]);
        y1 = Math.max(y1, l.points[2 * pt + 1]);
      }
      if (n) out.push({ x0, y0, x1, y1 });
    },
    glyphs() {},
    text() {},
    image(l: ImageList, i: number) {
      out.push({ x0: l.x[i], y0: l.y[i], x1: l.x[i] + l.w[i], y1: l.y[i] + l.h[i] });
    },
    clip() {},
    unclip() {},
  });
  return out;
}
