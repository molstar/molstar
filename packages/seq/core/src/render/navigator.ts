/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { type Issue } from '@molstar/core/repr/tree';
import { projectEnd, projectPoint } from '@molstar/seqx-core/resolve';
import { blockPositionAt, blockSpan, nearestOffsetOf, offsetOf, placesInto, positionAt } from '@molstar/seqx-core/view';
import type { AxisTransform } from './axis.js';
import { type RenderConfig } from './config.js';
import { type DrawLists, NoColor, type PassName, type PassWriter } from './draw-lists.js';
import { partOn } from './marks.js';
import type { Emphasis } from './layer-type.js';
import { BlockPlacement } from './layout/blocks.js';
import type { WindowPart } from './layout/types.js';
import type { SeqTokens } from './tokens.js';
import { hitView, type SectionLayout } from './tracks.js';
import { drawNavigatorWindow, floorExtent, NavWindow, NavWindowFlags, navigatorPart } from './window.js';

import {
  type NavigatorHit,
  blockIndexAt,
  blockUnitAtPosition,
  NavigatorWindow,
  partOf,
  windowAtFit,
  wantsStrip,
  rowsLayOf,
  ticks,
} from './navigator-parts.js';

export { blockPlaceText, type NavigatorHit, NavigatorWindow, navigatorBox, windowText } from './navigator-parts.js';

import { NavigatorBase } from './navigator-base.js';
import type { NavigatorFrame, NavigatorFrameInput } from './navigator-base.js';

export type { NavigatorFrame, NavigatorFrameInput, NavigatorRendererOptions } from './navigator-base.js';

// scratch and constants of the drawing, module-private (a draw reads no property of an imported object)
const ZoomHysteresisPx = 8;
const PlaceBlocks = new Int32Array(64);
const PlaceUnits = new Float64Array(64);
const FlagNoStart: number = NavWindowFlags.noStart;
const FlagNoEnd: number = NavWindowFlags.noEnd;
const BlockXs = new Float64Array(4);
const Ext = new Float64Array(8);
const MarkUs = new Float64Array(4);
const Line = new Float64Array(4);
const FlagMark: number = NavWindowFlags.outline | NavWindowFlags.quiet;
const None: number = NoColor;
const Seen: number[] = [];
const OwnPasses: readonly PassName[] = ['window', 'highlight', 'selection', 'hover', 'cursor'];

/**
 * Draws a navigator: its picture through a renderer of its own (the layout interface at the fit window), its windows
 * from the view's last layout. Allocates nothing once warm on a window move (the picture's key holds).
 */
export class NavigatorRenderer extends NavigatorBase {
  /** Lays the picture out (or keeps it), then draws the windows and marks over it. Valid until the next call. */
  frame(input: NavigatorFrameInput): NavigatorFrame {
    const f = this.frameRecord;
    // the bands and the rows fill the navigator's height: the picture is made again for another
    if (this.fills && input.size.height !== this.builtHeight) this.build(input.size.height);
    const r = this.renderer;
    if (!r) {
      f.rebuilt = false;
      f.base = false;
      f.wants = [];
      return f;
    }
    const fi = this.frameInput;
    this.boxHeight = input.size.height;
    fi.size.width = input.size.width;
    fi.size.height = input.size.height + this.zoomHeight;
    fi.dpr = input.dpr ?? 1;
    fi.settled = input.settled ?? false;
    fi.cvdSafe = input.cvdSafe ?? false;
    if (input.overlay) fi.overlay = input.overlay;
    else delete fi.overlay;
    this.setViewports(input.view);
    let frame = r.frame(fi);
    // the rows fit the band: their px a row follows the count the picture placed (a second layout, on a change only)
    if (frame.rebuilt && this.fitRows(frame.layout, input.view)) frame = r.frame(fi);
    if (frame.rebuilt && this.leaveOut(frame.layout)) return this.frame(input);
    this.layout = frame.layout;
    if (frame.rebuilt) {
      this.wants = frame.wants.filter((w) => w.priority !== 'margin');
      for (const i of frame.layout.issues) if (!this.issues.includes(i)) this.issues.push(i);
    }
    this.placeWindows(input.view);
    this.readReach();
    // the zoom strip (_spec/README.md "Render", Navigators; made again when it comes or goes) at the box's bottom
    // whatever the context's content; its viewport follows the span (laid out again)
    if (this.zoomFlips() || this.gapMoved()) {
      this.build();
      return this.frame(input);
    }
    if (
      this.zoomed &&
      (this.reach[0] !== this.zoomU0 || this.reach[1] !== this.zoomU1) &&
      this.reach[1] > this.reach[0]
    ) {
      this.zoomU0 = this.reach[0];
      this.zoomU1 = this.reach[1];
      this.setViewports(input.view, true);
      frame = r.frame(fi);
      this.layout = frame.layout;
      if (frame.rebuilt) this.wants = frame.wants.filter((w) => w.priority !== 'margin');
      this.placeWindows(input.view);
    }
    this.drawOver(frame.lists, input);
    f.lists = frame.lists;
    f.rebuilt = frame.rebuilt;
    f.base = frame.base;
    f.wants = this.wants;
    return f;
  }

  /**
   * Whether the zoom strip comes or goes: a positions navigator whose reachable span is narrower than `minPx` on the
   * picture (the context) has one. Its viewport is the span then.
   */
  private zoomFlips(): boolean {
    const ctx = this.layout?.sections[0]?.axis;
    const u0 = this.reach[0];
    const u1 = this.reach[1];
    const nav = this.navigator;
    // a span not known yet (the section's length pending) keeps what is drawn
    if (!ctx || !(u1 > u0)) return false;
    // a strip along the view only (a side part's length is the view's height, not its axis), a height to draw in
    const can = nav.axes === 'positions' && nav.place !== 'left' && nav.place !== 'right' && this.zoomTrack > 0;
    const span = (u1 - u0) * ctx.scale;
    const want = can && wantsStrip(this.overBlocks, this.zoomStrip, this.windowPx, span, this.zoomed, ZoomHysteresisPx);
    if (want === this.zoomed) return false;
    this.zoomed = want;
    this.zoomU0 = u0;
    this.zoomU1 = u1;
    for (const k of Object.keys(this.viewports)) delete this.viewports[k];
    this.viewportsFor = {};
    delete this.frameInput.viewports;
    return true;
  }

  /**
   * Whether the gap above the zoom strip moves: the strip sits at the box's bottom (the gap after it), whatever the
   * context's content takes (a strip's small tracks, its axis); a context taller than the box keeps the gap's least.
   */
  private gapMoved(): boolean {
    const ctx = this.layout?.sections[0];
    if (!this.zoomed || !ctx) return false;
    const gap = Math.max(0, this.boxHeight + this.zoomGap - (ctx.y + ctx.height));
    if (Math.abs(gap - this.sectionGap) < 1e-6) return false;
    this.sectionGap = gap;
    return true;
  }

  /** Over blocks: the picture's block under x, its system and the position there; undefined off the picture. */
  private blockHit(x: number): { block: number; system: string; position: number } | undefined {
    const work = this.work;
    const blocks = this.picture?.section.blocks;
    if (!work || !blocks?.length) return undefined;
    const p = work.placement;
    const u = this.unitOn(work, x);
    const k = p instanceof BlockPlacement && !p.stacked ? p.row.blockAt(x) : -1;
    const b = k >= 0 ? k : blockIndexAt(blocks, u, 'start');
    const position = blockPositionAt(blocks[b], u, 'start');
    return position === undefined ? undefined : { block: b, system: blocks[b].axis.system, position };
  }

  /**
   * Over blocks: the bound section's units of the picture's block under x (what a double-click fits, a switch shows
   * alone), its system and index; `start` and `end` NaN where the section shows none of it; undefined off a block.
   */
  blockAt(
    section: string,
    x: number,
    context = false,
  ): { start: number; end: number; system: string; block: number } | undefined {
    const blocks = this.picture?.section.blocks;
    const s = this.options.view.sections.find((v) => v.address === section);
    const work = context ? this.layout?.sections[0] : this.work;
    if (!this.overBlocks || !blocks || !s || !work) return undefined;
    const p = work.placement;
    const k = p instanceof BlockPlacement && !p.stacked ? p.row.blockAt(x) : -1;
    if (k < 0) return undefined;
    const block = blocks[k];
    const o = block.offset ?? 0;
    const a = this.fromPicture(s, o, 'start');
    const b = this.fromPicture(s, o + blockSpan(block), 'end');
    const shown = a !== undefined && b !== undefined && b > a;
    return { start: shown ? a : Number.NaN, end: shown ? b : Number.NaN, system: block.axis.system, block: k };
  }

  /** The windows (the top one last), then the marks, in the navigator's own passes. */
  private drawOver(lists: DrawLists, input: NavigatorFrameInput): void {
    lists.reset(OwnPasses);
    const out = lists.writer('window');
    const tokens = this.renderer!.tokens;
    const c = this.config.navigator.window;
    const dpr = input.dpr ?? 1;
    const moves = input.moves;
    if (this.zoomed) this.drawContext(out, dpr, c, tokens);
    // the reachable span's wash: fixed, where a window cannot go
    for (let i = 0; i < this.windows.length; i++) {
      const w = this.windows[i];
      if (!w.shown || !(w.axes & NavWindowFlags.x)) continue;
      const box = this.work!.placement.box;
      const wash = tokens.colors['seq-wash'];
      if (w.lo > box.x) {
        const g = out.args;
        g[0] = box.x;
        g[1] = w.ly;
        g[2] = w.lo - box.x;
        g[3] = w.hy - w.ly;
        g[4] = wash;
        g[5] = 0;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
      if (w.hi < box.x + box.width) {
        const g = out.args;
        g[0] = w.hi;
        g[1] = w.ly;
        g[2] = box.x + box.width - w.hi;
        g[3] = w.hy - w.ly;
        g[4] = wash;
        g[5] = 0;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        out.rectArgs();
      }
    }
    for (let pass = 0; pass < 2; pass++)
      for (let i = 0; i < this.windows.length; i++) {
        const w = this.windows[i];
        const isTop = w.section === input.top;
        if (!w.shown || (pass === 0) === isTop) continue;
        const m = moves ? (moves[w.section] ?? 0) : -1;
        // at the fit the frame is quiet like an idle one, no target: the row is a sweep (`windowAtFit`)
        const flags = m < 0 ? NavWindowFlags.outline : m === 0 || windowAtFit(w) ? NavWindowFlags.idle : m & w.axes;
        const lit = input.lit?.section === w.section ? input.lit.part : undefined;
        this.drawWindow(out, w, flags, m > 0 ? lit : undefined, dpr, c, tokens);
      }
    // a sweep: the region it would zoom to, an outline over the picture
    const sw = input.sweep;
    if (sw) {
      const W = NavWindow;
      W[0] = W[4] = Math.min(sw.x0, sw.x1);
      W[1] = W[5] = Math.max(sw.x0, sw.x1);
      const work = this.work!;
      W[2] = W[6] = W[10] = this.zoomed ? work.y : 0;
      W[3] = W[7] = W[11] = this.zoomed ? work.y + work.height : this.frameInput.size.height;
      W[8] = W[0];
      W[9] = W[1];
      drawNavigatorWindow(out, NavWindowFlags.outline, 'body', dpr, c, tokens);
    }
    this.drawMarks(lists, input, tokens, dpr);
  }

  /**
   * The context of a zoom strip: the reachable span's mark on the picture (an outline at least `window.minPx` wide,
   * centred on the truth, its hairline at the truth; no grips, no wash) and two connector hairlines from its ends to the
   * zoom strip's ends.
   */
  private drawContext(out: PassWriter, dpr: number, c: RenderConfig['navigator']['window'], tokens: SeqTokens): void {
    const ctx = this.layout!.sections[0];
    const work = this.work!;
    const axis = ctx.axis;
    const m = this.markX;
    m[0] = m[1] = m[2] = m[3] = Number.NaN;
    if (!axis || !(this.reach[1] > this.reach[0])) return;
    const t = axis.total;
    // over a cyclic origin two pieces: the arc's start to the axis's end, its start to the arc's end
    // inline, never `this.wraps`: its `turn` getter returns a genome's total (past a small integer: boxed, 16 B a frame)
    const wraps = t > 0 && this.reach[1] > t + 1e-9;
    // the ends through `MarkUs` (a number chosen by a conditional is boxed when passed on)
    const U = MarkUs;
    U[0] = this.reach[0];
    U[1] = wraps ? t : this.reach[1];
    U[2] = 0;
    U[3] = this.reach[1] - t;
    for (let k = 0; k < (wraps ? 2 : 1); k++) {
      this.xIntoFrom(ctx, U, 2 * k, BlockXs, 0);
      this.xIntoFrom(ctx, U, 2 * k + 1, BlockXs, 1);
      // the piece's ends through `BlockXs` (two fractional numbers passed to a call that is not inlined are boxed)
      this.drawMark(out, ctx, k, dpr, c, tokens);
    }
    // the connectors: from the mark's foot to the zoom strip's top, its start to the strip's start, its end to its end
    // (over the origin the first piece's start, the second's end)
    if (this.config.navigator.zoomStrip.funnel === 'none') return;
    const g = out.args;
    const wb = work.placement.box;
    for (let k = 0; k < 2; k++) {
      g[4] = None;
      g[5] = ((tokens.colors['seq-window'] & 0xffffff00) | Math.round(c.idleAlpha * 255)) >>> 0;
      g[6] = c.strokePx;
      g[7] = 0;
      // display units run left to right on both: the arc's start is the strip's start
      Line[0] = k === 0 ? m[0] : wraps ? m[3] : m[1];
      Line[1] = ctx.y + ctx.height;
      Line[2] = k === 0 ? wb.x : wb.x + wb.width;
      Line[3] = work.y;
      out.pathArgs(Line, 2, false);
    }
  }

  /**
   * A piece of the mark between `BlockXs[0]` and `BlockXs[1]` on the context: an outline at least `window.minPx` wide, always centred on
   * the truth (whatever `window.floor`: the mark says where) and clamped inside the picture, so the truth is inside it;
   * its hairline at the truth. Its drawn ends into `markX` (piece `k`).
   */
  private drawMark(
    out: PassWriter,
    ctx: SectionLayout,
    k: number,
    dpr: number,
    c: RenderConfig['navigator']['window'],
    tokens: SeqTokens,
  ): void {
    const box = ctx.placement.box;
    const a = BlockXs[0];
    const b = BlockXs[1];
    Ext[0] = a < b ? a : b;
    Ext[1] = a < b ? b : a;
    Ext[2] = box.x;
    Ext[3] = box.x + box.width;
    floorExtent(Ext, 0, this.markConfig(c));
    this.markX[2 * k] = Ext[0];
    this.markX[2 * k + 1] = Ext[1];
    const W = NavWindow;
    W[0] = W[8] = Ext[0];
    W[1] = W[9] = Ext[1];
    W[4] = Math.min(a, b);
    W[5] = Math.max(a, b);
    W[2] = W[6] = W[10] = ctx.y;
    W[3] = W[7] = W[11] = ctx.y + ctx.height;
    drawNavigatorWindow(out, FlagMark, undefined, dpr, c, tokens);
  }

  /** `window` with its floor `centred` (the mark's), made once. */
  private centred: RenderConfig['navigator']['window'] | undefined;
  private markConfig(c: RenderConfig['navigator']['window']): RenderConfig['navigator']['window'] {
    if (c.floor === 'centred') return c;
    if (this.centred?.minPx !== c.minPx) this.centred = { ...c, floor: 'centred' };
    return this.centred;
  }

  private drawWindow(
    out: PassWriter,
    w: NavigatorWindow,
    flags: number,
    lit: WindowPart | undefined,
    dpr: number,
    c: RenderConfig['navigator']['window'],
    tokens: SeqTokens,
  ): void {
    const W = NavWindow;
    W[2] = w.dy0;
    W[3] = w.dy1;
    W[6] = w.y0;
    W[7] = w.y1;
    W[8] = w.lo;
    W[9] = w.hi;
    W[10] = w.ly;
    W[11] = w.hy;
    if (w.wrap) {
      // over the origin: its two pieces, the wash between them only, a grip at each true end (the first piece's start,
      // the second's end): where a piece meets the axis's end it is body
      W[0] = w.w0;
      W[1] = w.w1;
      W[4] = w.w0;
      W[5] = w.w1;
      W[8] = w.w0;
      W[9] = w.w1;
      drawNavigatorWindow(out, flags | FlagNoStart, lit, dpr, c, tokens);
      W[8] = w.w1;
      W[9] = w.x0;
    }
    W[0] = w.wrap ? w.x0 : w.dx0;
    W[1] = w.wrap ? w.x1 : w.dx1;
    W[4] = w.x0;
    W[5] = w.x1;
    drawNavigatorWindow(out, w.wrap ? flags | FlagNoEnd : flags, lit, dpr, c, tokens);
  }

  /** The selection's and the highlights' ticks at their places, the cursor's hairline, the pointer's hairline. */
  private drawMarks(lists: DrawLists, input: NavigatorFrameInput, tokens: SeqTokens, dpr: number): void {
    const pic = this.work;
    const axis = pic?.axis;
    if (!pic || !axis) return;
    const marks = this.config.navigator.marks;
    const e = input.emphasis;
    const system = axis.system;
    const box = pic.placement.box;
    const h = Math.min(4, box.height);
    const y = box.y + box.height - h;
    const px = 1 / dpr;
    if (this.overBlocks) this.drawBlockMarks(lists, pic, e, marks, y, h, box, tokens, px);
    else if (e && marks.includes('selection')) {
      const part = partOn(e.selection, system) ?? partOn(e.projected?.selection, system);
      if (part) ticks(lists.writer('selection'), axis, part.ranges, y, h, px, tokens.colors['msui-accent']);
    }
    if (!this.overBlocks && e?.highlights && marks.includes('highlights')) {
      const out = lists.writer('highlight');
      for (let k = 0; k < e.highlights.length; k++) {
        const part = partOn(e.highlights[k].location, system);
        if (part) ticks(out, axis, part.ranges, y, h, px, e.highlights[k].color);
      }
    }
    if (!this.overBlocks && e?.cursor && e.cursor.system === system && marks.includes('cursor')) {
      const x = axis.x(e.cursor.position);
      if (Number.isFinite(x)) {
        const g = lists.writer('cursor').args;
        g[0] = x - px / 2;
        g[1] = box.y;
        g[2] = px;
        g[3] = box.height;
        g[4] = tokens.colors['msui-text'];
        g[5] = 0;
        g[6] = 0;
        g[7] = 0;
        g[8] = 0;
        lists.writer('cursor').rectArgs();
      }
    }
    // the pointer's hairlines: down the picture at its position (none on a rows navigator: it has no positions axis;
    // with a zoom strip, down the part it is over), across it at its row
    const px0 = input.pointerX ?? Number.NaN;
    if (this.navigator.axes !== 'rows' && Number.isFinite(px0) && px0 >= box.x && px0 < box.x + box.width) {
      const over = this.zoomed && (input.pointerY ?? Number.NaN) < pic.y ? this.layout!.sections[0] : pic;
      const top = this.zoomed ? over.y : box.y;
      const bottom = this.zoomed ? over.y + over.height : this.frameInput.size.height;
      lists.writer('hover').fill(px0 - px / 2, top, px, bottom - top, tokens.colors['seq-hover']);
    }
    const py0 = input.pointerY ?? Number.NaN;
    if (Number.isFinite(py0) && this.navigator.axes !== 'positions') {
      const g = lists.writer('hover').args;
      g[0] = box.x;
      g[1] = py0 - px / 2;
      g[2] = box.width;
      g[3] = px;
      g[4] = tokens.colors['seq-hover'];
      g[5] = 0;
      g[6] = 0;
      g[7] = 0;
      g[8] = 0;
      lists.writer('hover').rectArgs();
    }
  }

  /**
   * Over blocks: the selection's and the highlights' ticks on each block showing them, the keyboard cursor's hairline at
   * every place of its position (overlapping blocks: each).
   */
  private drawBlockMarks(
    lists: DrawLists,
    pic: SectionLayout,
    e: Emphasis | undefined,
    marks: RenderConfig['navigator']['marks'],
    y: number,
    h: number,
    box: { readonly y: number; readonly height: number },
    tokens: SeqTokens,
    px: number,
  ): void {
    const blocks = this.picture?.section.blocks;
    if (!e || !blocks) return;
    const tick = (out: PassWriter, loc: Emphasis['selection'], color: number) => {
      for (let b = 0; b < blocks.length; b++) {
        const part = partOn(loc, blocks[b].axis.system);
        if (!part) continue;
        const o = blocks[b].offset ?? 0;
        for (let k = 0; k + 1 < part.ranges.length; k += 2) {
          const ua = blockUnitAtPosition(blocks[b], part.ranges[k], o);
          const ub = blockUnitAtPosition(blocks[b], part.ranges[k + 1], o);
          if (ua === undefined || ub === undefined) continue;
          const xa = this.xOn(pic, ua);
          const xb = this.xOn(pic, ub);
          out.fill(Math.min(xa, xb), y, Math.max(px, Math.abs(xb - xa)), h, color);
        }
      }
    };
    if (marks.includes('selection')) {
      const sel = e.selection ?? e.projected?.selection;
      if (sel) tick(lists.writer('selection'), sel, tokens.colors['msui-accent']);
    }
    if (e.highlights && marks.includes('highlights'))
      for (let k = 0; k < e.highlights.length; k++)
        tick(lists.writer('highlight'), e.highlights[k].location, e.highlights[k].color);
    if (e.cursor && marks.includes('cursor')) {
      const n = placesInto(this.picture!.section, e.cursor.system, e.cursor.position, PlaceBlocks, PlaceUnits);
      for (let k = 0; k < n; k++) {
        const x = this.xOn(pic, PlaceUnits[k]);
        if (Number.isFinite(x))
          lists.writer('cursor').fill(x - px / 2, box.y, px, box.height, tokens.colors['msui-text']);
      }
    }
  }

  /** The part under a point of the navigator, the position, row and band there; undefined off the picture. */
  hit(x: number, y: number, touch = false, top?: string): NavigatorHit | undefined {
    const pic = this.work;
    if (!pic) return undefined;
    const c = this.config.navigator.window;
    const nav = this.navigator;
    // on the context, the zoom strip's: no target; over the mark it names the region. The gap below it is the strip's
    // (a finger a few px above the strip lands on it)
    const ctx = this.layout!.sections[0];
    if (this.zoomed && y < ctx.y + ctx.height) return this.contextHit(x, y);
    let best: { w: NavigatorWindow; part: WindowPart; axis: 'positions' | 'rows' } | undefined;
    Seen.length = 0;
    for (let k = -1; k < this.windows.length; k++) {
      // the window the pointer last acted on first (overlapping windows: it takes the press)
      const i = k < 0 ? this.windows.findIndex((w) => w.section === top) : k;
      if (i < 0 || Seen.includes(i)) continue;
      Seen.push(i);
      const w = this.windows[i];
      if (!w.shown) continue;
      const along = w.axes & NavWindowFlags.x ? partOf(x, w, c, touch) : 'body';
      const across = w.axes & NavWindowFlags.y ? navigatorPart(y, w.dy0, w.dy1, w.ly, w.hy, c, touch) : 'body';
      if (along === undefined && across === undefined) continue;
      const a = along ?? 'outside';
      const b = across ?? 'outside';
      let part: WindowPart;
      let axis: 'positions' | 'rows' = nav.axes === 'rows' ? 'rows' : 'positions';
      if ((b === 'start' || b === 'end') && a !== 'outside') {
        part = b;
        axis = 'rows';
      } else if ((a === 'start' || a === 'end') && b !== 'outside') part = a;
      else if (a === 'body' && b === 'body') part = 'body';
      else part = 'outside';
      if (part !== 'outside') {
        best = { w, part, axis };
        break;
      }
      best ??= { w, part, axis };
    }
    const section = best?.w.section ?? nav.of[0];
    if (section === undefined) return undefined;
    const h = this.layout ? hitView(this.layout, x, y) : undefined;
    const axis = pic.axis;
    // over blocks: the block under the point, its system and the position there
    const bh = this.overBlocks ? this.blockHit(x) : undefined;
    const p = bh
      ? bh.position
      : axis && !this.overBlocks
        ? axis.positionAt(Math.max(axis.x0, Math.min(axis.x0 + axis.width - 1e-9, x)))
        : -1;
    const position = p >= 0 ? p : undefined;
    const row = this.rowAt(y);
    return {
      section,
      part: best?.part ?? 'outside',
      axis: best?.axis ?? (nav.axes === 'rows' ? 'rows' : 'positions'),
      ...(position !== undefined && nav.axes !== 'rows' ? { position } : {}),
      ...(bh ? { system: bh.system, block: bh.block } : {}),
      ...(row !== undefined ? { row } : {}),
      ...(h?.state !== undefined ? { band: h.state } : {}),
    };
  }

  /** A point on the context of a zoom strip: its position and band, over the mark or not; never a window's part. */
  private contextHit(x: number, y: number): NavigatorHit | undefined {
    const section = this.navigator.of[0];
    const ctx = this.layout?.sections[0];
    if (section === undefined || !ctx) return undefined;
    const h = hitView(this.layout!, x, y);
    const axis = ctx.axis;
    const blocks = this.overBlocks ? this.picture?.section.blocks : undefined;
    let p = axis && !blocks ? axis.positionAt(Math.max(axis.x0, Math.min(axis.x0 + axis.width - 1e-9, x))) : -1;
    let at: { system: string; block: number } | undefined;
    if (blocks?.length) {
      const u = this.unitOn(ctx, x);
      const b = blockIndexAt(blocks, u, 'start');
      p = blockPositionAt(blocks[b], u, 'start') ?? -1;
      at = { system: blocks[b].axis.system, block: b };
    }
    return {
      section,
      part: 'outside',
      axis: 'positions',
      context: this.onMark(x) ? 'mark' : 'outside',
      ...(p >= 0 ? { position: p } : {}),
      ...(at ?? {}),
      ...(h?.state !== undefined ? { band: h.state } : {}),
    };
  }

  /** Whether x is over a piece of the mark as last drawn. */
  private onMark(x: number): boolean {
    const m = this.markX;
    return (x >= m[0] && x <= m[1]) || (x >= m[2] && x <= m[3]);
  }

  /** The display row (below the pinned band, fractional) at y on the picture's rows; undefined with none. */
  rowAt(y: number): number | undefined {
    const nr = this.layout ? rowsLayOf(this.layout.sections[0]) : undefined;
    if (!nr || nr.total === undefined || !(nr.rowScale > 0)) return undefined;
    const top = nr.box.y + nr.pinnedHeight;
    const rows = nr.display.count - nr.display.pinned;
    return Math.max(0, Math.min(rows, (y - top) / nr.rowScale));
  }

  /**
   * The bound section's display unit at x (its own axis's offset), through the inverse of the path from its system:
   * what a drag on its window moves it to. Past the picture's ends the drag keeps going at the navigator's scale, as the
   * path carries the edge (a codon a protein position). Undefined while the picture is not placed or the place maps
   * nowhere.
   */
  sectionUnitAt(section: string, x: number, context = false): number | undefined {
    const axis = this.work?.axis;
    if (!axis) return undefined;
    if (this.overBlocks) {
      // over blocks the context is the genome: a click there is a jump (`context`: read on it, not the strip)
      const on = context ? this.layout!.sections[0] : this.work!;
      const s = this.options.view.sections.find((v) => v.address === section);
      return s ? this.fromPicture(s, this.unitOn(on, x), 'start') : undefined;
    }
    const u = axis.u0 + (x - axis.x0) / axis.scale;
    const inside = Math.max(0, Math.min(axis.total, u));
    const o = this.unitOfPicture(section, axis, inside);
    if (o === undefined || inside === u) return o;
    const step = u < 0 ? 1 : -1;
    const near = this.unitOfPicture(section, axis, inside + step);
    return near === undefined ? o : o + (u - inside) * ((near - o) / step);
  }

  /** The bound section's display unit at the picture's unit `u` (within the picture). */
  private unitOfPicture(section: string, axis: AxisTransform, u: number): number | undefined {
    const s = this.options.view.sections.find((v) => v.address === section);
    const path = this.paths.get(section);
    if (!s || !path) return undefined;
    // the picture's far end is a window's end edge: carried as one
    const end = u >= axis.total;
    let p = positionAt(axis, u, end ? 'end' : 'start');
    if (p === undefined) return undefined;
    if (path.back) {
      p = end ? projectEnd(path.back, p) : projectPoint(path.back, p);
      if (p === undefined) return undefined;
    }
    return offsetOf(s.axis, p) ?? nearestOffsetOf(s.axis, p);
  }

  /** x on the navigator of a bound section's display unit (through the path from its system): a sweep's ends. */
  xOfSectionUnit(section: string, u: number): number {
    const axis = this.work?.axis;
    const s = this.options.view.sections.find((v) => v.address === section);
    const path = this.paths.get(section);
    if (!axis || !s || !path) return Number.NaN;
    if (this.overBlocks) {
      const pu = this.toPicture(s, u, 'start');
      return pu === undefined ? Number.NaN : this.xOn(this.work!, pu);
    }
    let p = positionAt(s.axis, u);
    if (p === undefined) return Number.NaN;
    if (path.to) p = projectPoint(path.to, p);
    if (p === undefined) return Number.NaN;
    const o = offsetOf(axis, p) ?? nearestOffsetOf(axis, p);
    return o === undefined ? Number.NaN : axis.xOfUnit(this.onStrip(o));
  }

  /** The issues the navigator met: its bindings, its layers left out, its picture's. */
  get allIssues(): readonly Issue[] {
    return this.issues;
  }
}
