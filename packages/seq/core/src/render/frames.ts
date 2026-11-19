/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// The `frame` layer type (_spec/README.md "Render", Layer types): an annotation's regions outlined across tracks of its
// section (ESPript's frames, EC-100). Strip scope: on each strip a region's outline covers `ctx.span` (the bands of the
// tracks `across` names, its own when it names none) over the region's x; an end the strip cuts (a line break, the
// view's edge, a piece's) is left open, as every cut end is drawn square (the frame rule); its title (the item's
// label) is drawn in the layer's own band where the region really starts. A strip with none of its regions leaves its
// band unused (text collapses it there). Its outline and title hit as the region's item.

import type { ResolvedLayer } from '@molstar/seqx-core/view';
import type { StaticModel } from '@molstar/seqx-core/resolve';
import { Cut } from './axis.js';
import { Anchor, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { beginItems, type ItemsBase, type ItemsPrep, itemHit, newItemsBase, readItems } from './items.js';
import { layerSource, toSection } from './layer-data.js';
import type { Box, Hit, LayerType, LayoutContext } from './layer-type.js';

// the anchor a number read once: a draw reads no property of an imported object (_spec/README.md "Render", Draw lists)
const AnchorStart: number = Anchor.start;
const CutLeft: number = Cut.left;
const CutRight: number = Cut.right;

/** How near an outline a point hits it, CSS px (fixed, as other hit boxes are). */
const Slop = 3;

export interface FramePrep extends ItemsPrep {}

export interface FrameLay extends ItemsBase {
  /** The outlines on this strip: `x0, x1, cut, row` (the row in `rows`), and whether its title is drawn here. */
  m: number;
  marks: number[];
  titled: Uint8Array;
  /** What the outlines cover (the named tracks' bands on the strip) and the title's baseline. */
  span: Box;
  font: string;
  titleY: number;
}

const Pieces: number[] = [];
const Spans: number[] = [];

export const FrameType: LayerType<FramePrep, FrameLay> = {
  name: 'frame',

  prepare(layer: ResolvedLayer, model: StaticModel): FramePrep {
    const src = layerSource(model, layer);
    return { ...(src ? { src } : {}), owner: layer.data?.address ?? '' };
  },

  layout(prep, ctx: LayoutContext, data, _overlay, prev?: FrameLay): FrameLay {
    const lay: FrameLay = prev ?? {
      ...newItemsBase(ctx, prep),
      m: 0,
      marks: [],
      titled: new Uint8Array(8),
      span: ctx.box,
      font: ctx.config.itemFont,
      titleY: 0,
    };
    beginItems(lay, ctx, prep);
    lay.m = 0;
    lay.marks.length = 0;
    lay.span = ctx.span && ctx.span.height > 0 ? ctx.span : ctx.box;
    lay.font = ctx.config.itemFont;
    lay.titleY = ctx.box.y + ctx.box.height - 2;
    const src = prep.src;
    if (!src) return lay;
    readItems(lay, ctx, data, src, ['label']);
    const axis = ctx.axis;
    const rows = lay.rows;
    for (let i = 0; i < rows.n; i++) {
      toSection(src, rows.starts[i], rows.ends[i], Pieces);
      for (let k = 0; k + 1 < Pieces.length; k += 2) {
        Spans.length = 0;
        axis.spans(Pieces[k], Pieces[k + 1], Spans);
        for (let j = 0; j + 2 < Spans.length; j += 3) {
          lay.marks.push(Spans[j], Spans[j + 1], Spans[j + 2], i);
          lay.m++;
        }
      }
    }
    if (lay.titled.length < lay.m) lay.titled = new Uint8Array(lay.m * 2);
    // a title where the region starts: its first mark in display order whose left end is its own
    const seen = new Set<number>();
    for (let j = 0; j < lay.m; j++) {
      const row = lay.marks[4 * j + 3];
      lay.titled[j] = 0;
      if ((lay.marks[4 * j + 2] & CutLeft) === 0 && !seen.has(row)) {
        seen.add(row);
        lay.titled[j] = 1;
      }
    }
    return lay;
  },

  // its band holds the titles; a strip with an outline keeps it (a band collapsed is not laid out, so its outlines
  // would go with it)
  used(lay) {
    return lay.m > 0 ? lay.box.height : 0;
  },

  draw(lay, pass: PassName, out: PassWriter) {
    if (pass !== 'data' && pass !== 'labels') return;
    const ink = lay.tokens.colors['msui-text-muted'];
    const g = out.args;
    const px = 1 / lay.axis.dpr;
    const y0 = lay.span.y;
    const y1 = lay.span.y + lay.span.height;
    for (let j = 0; j < lay.m; j++) {
      const x0 = lay.marks[4 * j];
      const x1 = lay.marks[4 * j + 1];
      const cut = lay.marks[4 * j + 2];
      if (pass === 'labels') {
        if (!lay.titled[j]) continue;
        const text = lay.rows.labels[lay.marks[4 * j + 3]] ?? String(lay.rows.keys[lay.marks[4 * j + 3]]);
        g[0] = x0;
        g[1] = lay.titleY;
        g[2] = 0;
        g[3] = Math.max(0, x1 - x0);
        g[4] = ink;
        out.textArgs(text, AnchorStart, lay.font);
        continue;
      }
      // the top and the bottom always; a side only where the region really ends (a cut end is open)
      line(out, x0, y0, x1 - x0, px, ink);
      line(out, x0, y1 - px, x1 - x0, px, ink);
      if ((cut & CutLeft) === 0) line(out, x0, y0, px, y1 - y0, ink);
      if ((cut & CutRight) === 0) line(out, x1 - px, y0, px, y1 - y0, ink);
    }
  },

  hit(lay, x, y): Hit | undefined {
    const y0 = lay.span.y;
    const y1 = lay.span.y + lay.span.height;
    for (let j = lay.m - 1; j >= 0; j--) {
      const x0 = lay.marks[4 * j];
      const x1 = lay.marks[4 * j + 1];
      const cut = lay.marks[4 * j + 2];
      const inX = x >= x0 - Slop && x <= x1 + Slop;
      const onTitle = lay.titled[j] === 1 && y >= lay.box.y && y < lay.box.y + lay.box.height && x >= x0 && x < x1;
      const left = (cut & CutLeft) === 0 && Math.abs(x - x0) <= Slop;
      const right = (cut & CutRight) === 0 && Math.abs(x - x1) <= Slop;
      const onSide = y >= y0 - Slop && y <= y1 + Slop && (left || right);
      const onEdge = inX && (Math.abs(y - y0) <= Slop || Math.abs(y - y1) <= Slop);
      if (!onTitle && !onSide && !onEdge) continue;
      const row = lay.marks[4 * j + 3];
      const base = { layer: lay.key, section: lay.section, system: lay.system };
      return itemHit(base, lay, lay.rows.keys[row], lay.rows.starts[row], lay.rows.ends[row]);
    }
    return undefined;
  },
};

/** A hairline as a filled rect, through `args` (nothing boxed). */
function line(out: PassWriter, x: number, y: number, w: number, h: number, ink: number) {
  const g = out.args;
  g[0] = x;
  g[1] = y;
  g[2] = w;
  g[3] = h;
  g[4] = ink;
  g[5] = NoColor;
  g[6] = 0;
  g[7] = 0;
  g[8] = 0;
  out.rectArgs();
}
