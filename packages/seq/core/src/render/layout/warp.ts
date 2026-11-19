/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Warps from a strip's frame to the page (_spec/README.md "Render", Section layouts): what a layer draws in its strip's
// frame reaches the page through its strip's warp, so no layer type knows a layout. A translation moves what a strip's
// draw added, in place (text: a line's lays kept across a scroll, only the translation changes); the polar warp takes
// what a strip drew into a scratch pass and writes it bent onto the page: a rect an annular sector (its squircle
// corners warped with it), a path's lines cut every few px so they bend, a glyph run's letters placed one by one with
// their angles, an image refused (its layout says so), a clip applied in the strip's frame before warping. Buffers are
// kept: nothing allocated once warm.

import { grown, NoColor, Op, PassList, Passes, PassWriter, PathCommand, type PathList } from '../draw-lists.js';
import { squirclePath } from '../squircle.js';

/** What a translation moves: the counts of a pass's lists when a strip began to draw. */
export class TranslateWarp {
  private rects = 0;
  private points = 0;
  private runs = 0;
  private letters = 0;
  private texts = 0;
  private images = 0;
  private clips = 0;

  /** Marks where a strip's draw begins in `out`'s pass. */
  open(out: PassWriter): void {
    const p = out.pass;
    this.rects = p.rects.n;
    this.points = p.paths.pointCount;
    this.runs = p.glyphs.n;
    this.letters = p.glyphs.letters;
    this.texts = p.texts.n;
    this.images = p.images.n;
    this.clips = p.clips.n;
  }

  /** Moves what was drawn since `open` by `dx`, `dy`. */
  close(out: PassWriter, dx: number, dy: number): void {
    if (dx === 0 && dy === 0) return;
    const p = out.pass;
    const r = p.rects;
    for (let i = this.rects; i < r.n; i++) {
      r.x[i] += dx;
      r.y[i] += dy;
    }
    const pts = p.paths.points;
    for (let i = this.points; i < p.paths.pointCount; i++) {
      pts[2 * i] += dx;
      pts[2 * i + 1] += dy;
    }
    const g = p.glyphs;
    for (let i = this.runs; i < g.n; i++) g.y[i] += dy;
    for (let j = this.letters; j < g.letters; j++) {
      g.x[j] += dx;
      g.ly[j] += dy;
    }
    const t = p.texts;
    for (let i = this.texts; i < t.n; i++) {
      t.x[i] += dx;
      t.y[i] += dy;
    }
    const im = p.images;
    for (let i = this.images; i < im.n; i++) {
      im.x[i] += dx;
      im.y[i] += dy;
    }
    const c = p.clips;
    for (let i = this.clips; i < c.n; i++) {
      c.x[i] += dx;
      c.y[i] += dy;
    }
  }
}

/**
 * The polar warp of a strip: x along the strip's frame to an angle (`a0` at `x0`, `sweep` radians over `width`,
 * clockwise from the top; negative counterclockwise), y down the frame to a radius (`r0` at y 0, inward). Letters are
 * upright or along the ring (`tangent`). A strip draws into `writer(out)` (a scratch pass sharing the lists' strings)
 * and `close(out)` writes it onto the page.
 */
export class PolarWarp {
  cx = 0;
  cy = 0;
  r0 = 0;
  x0 = 0;
  width = 1;
  a0 = 0;
  sweep = 2 * Math.PI;
  /** Letters along the ring (`tangent`) or upright. */
  tangent = false;
  /** The most px a cut segment spans along the ring. */
  step = 4;
  /** Primitives refused (an image): its layout reports them. */
  refused = 0;
  /** A scratch pass a pass, named as it (what a strip draws reads its pass's name: the emphasis does). */
  private readonly scratches = Passes.map((name) => new PassList(name));
  private readonly writers: (PassWriter | undefined)[] = [];
  private scratch = this.scratches[0];
  /** The clip rects open in the strip's frame (x, y, w, h each). */
  private clipStack = new Float64Array(16);
  private clipDepth = 0;
  /** The square's outline a rect is traced through. */
  private readonly square = new PassList('data');

  /** The page point of a frame point: into `px`, `py`. */
  px = 0;
  py = 0;
  map(x: number, y: number): void {
    const t = this.a0 + ((x - this.x0) / this.width) * this.sweep;
    const r = this.r0 - y;
    this.px = this.cx + r * Math.sin(t);
    this.py = this.cy - r * Math.cos(t);
  }

  /** The angle at a frame x. */
  angleAt(x: number): number {
    return this.a0 + ((x - this.x0) / this.width) * this.sweep;
  }

  /** The frame point of a page point (the inverse): into `px`, `py`; false at the center. */
  unmap(x: number, y: number): boolean {
    const dx = x - this.cx;
    const dy = y - this.cy;
    const r = Math.sqrt(dx * dx + dy * dy);
    if (r === 0) return false;
    let t = Math.atan2(dx, -dy) - this.a0;
    const full = 2 * Math.PI;
    if (this.sweep < 0) t = -t;
    t = ((t % full) + full) % full;
    this.px = this.x0 + (t / Math.abs(this.sweep)) * this.width;
    this.py = this.r0 - r;
    return true;
  }

  /** The writer a strip draws into for `out`'s pass. */
  writer(out: PassWriter): PassWriter {
    const k = out.pass.index;
    this.scratch = this.scratches[k];
    this.scratch.reset();
    this.clipDepth = 0;
    let w = this.writers[k];
    if (!w || w.lists !== out.lists) w = this.writers[k] = new PassWriter(out.lists, this.scratch);
    return w;
  }

  /** Writes what the strip drew onto `out`'s pass, warped. */
  close(out: PassWriter): void {
    const p = this.scratch;
    let r = 0;
    let pa = 0;
    let g = 0;
    let t = 0;
    let im = 0;
    let c = 0;
    for (let k = 0; k < p.opCount; k++) {
      switch (p.ops[k]) {
        case Op.rect:
          this.rect(out, p, r++);
          break;
        case Op.path:
          this.path(out, p.paths, pa++);
          break;
        case Op.glyphs:
          this.glyphs(out, p, g++);
          break;
        case Op.text:
          this.text(out, p, t++);
          break;
        case Op.image:
          im++;
          this.refused++;
          break;
        case Op.clip: {
          const l = p.clips;
          this.pushClip(l.x[c], l.y[c], l.w[c], l.h[c]);
          c++;
          break;
        }
        case Op.unclip:
          if (this.clipDepth > 0) this.clipDepth--;
          break;
      }
    }
    p.reset();
  }

  private pushClip(x: number, y: number, w: number, h: number) {
    let x0 = x;
    let y0 = y;
    let x1 = x + w;
    let y1 = y + h;
    if (this.clipDepth > 0) {
      const s = this.clipStack;
      const b = 4 * (this.clipDepth - 1);
      x0 = Math.max(x0, s[b]);
      y0 = Math.max(y0, s[b + 1]);
      x1 = Math.min(x1, s[b] + s[b + 2]);
      y1 = Math.min(y1, s[b + 1] + s[b + 3]);
    }
    this.clipStack = grown(this.clipStack, 4 * (this.clipDepth + 1));
    const b = 4 * this.clipDepth++;
    this.clipStack[b] = x0;
    this.clipStack[b + 1] = y0;
    this.clipStack[b + 2] = Math.max(0, x1 - x0);
    this.clipStack[b + 3] = Math.max(0, y1 - y0);
  }

  /** A rect: clipped in the frame, traced as its squircle outline, written as a bent path. */
  private rect(out: PassWriter, p: PassList, i: number) {
    const l = p.rects;
    let x = l.x[i];
    let y = l.y[i];
    let w = l.w[i];
    let h = l.h[i];
    if (this.clipDepth > 0) {
      const s = this.clipStack;
      const b = 4 * (this.clipDepth - 1);
      const x1 = Math.min(x + w, s[b] + s[b + 2]);
      const y1 = Math.min(y + h, s[b + 1] + s[b + 3]);
      x = Math.max(x, s[b]);
      y = Math.max(y, s[b + 1]);
      w = x1 - x;
      h = y1 - y;
      if (w <= 0 || h <= 0) return;
    }
    const q = this.square;
    q.reset();
    q.paths.begin();
    squirclePath(q.paths, x, y, w, h, l.radius[i], l.round[i]);
    q.paths.end(l.fill[i], l.stroke[i], l.strokeWidth[i], 0);
    this.path(out, q.paths, 0);
  }

  /** A path bent onto the page: its lines cut every `step` px of arc, its curves' points mapped. */
  private path(out: PassWriter, l: PathList, i: number) {
    out.pass.op(Op.path);
    const o = out.pass.paths;
    o.begin();
    let pt = l.pointStart[i];
    let lx = 0;
    let ly = 0;
    for (let c = l.cmdStart[i], end = c + l.cmdCount[i]; c < end; c++) {
      const cmd = l.commands[c];
      if (cmd === PathCommand.move) {
        lx = l.points[2 * pt];
        ly = l.points[2 * pt + 1];
        pt++;
        this.map(lx, ly);
        o.command(PathCommand.move);
        o.point(this.px, this.py);
      } else if (cmd === PathCommand.line) {
        const x = l.points[2 * pt];
        const y = l.points[2 * pt + 1];
        pt++;
        // the segment's length along the ring at its radius, cut every `step` px so it bends
        const r = Math.max(0, this.r0 - (ly + y) / 2);
        const arc = Math.abs(((x - lx) / this.width) * this.sweep) * r;
        const n = Math.max(1, Math.min(512, Math.ceil(arc / this.step)));
        for (let k = 1; k <= n; k++) {
          this.map(lx + ((x - lx) * k) / n, ly + ((y - ly) * k) / n);
          o.command(PathCommand.line);
          o.point(this.px, this.py);
        }
        lx = x;
        ly = y;
      } else if (cmd === PathCommand.cubic) {
        o.command(PathCommand.cubic);
        for (let k = 0; k < 3; k++, pt++) {
          this.map(l.points[2 * pt], l.points[2 * pt + 1]);
          o.point(this.px, this.py);
        }
        lx = l.points[2 * pt - 2];
        ly = l.points[2 * pt - 1];
      } else o.command(PathCommand.close);
    }
    o.end(l.fill[i], l.stroke[i], l.strokeWidth[i], l.dash[i]);
  }

  /** A glyph run placed letter by letter: each at its point, upright or along the ring. */
  private glyphs(out: PassWriter, p: PassList, i: number) {
    const l = p.glyphs;
    const placed = l.placed[i] === 1;
    out.pass.op(Op.glyphs);
    const o = out.pass.glyphs;
    const k = o.nextRun();
    o.font[k] = l.font[i];
    o.y[k] = l.y[i];
    o.color[k] = l.color[i];
    o.placed[k] = 1;
    const cs = this.clipDepth > 0 ? 4 * (this.clipDepth - 1) : -1;
    for (let j = l.start[i], end = j + l.count[i]; j < end; j++) {
      const x = l.x[j];
      if (cs >= 0 && (x < this.clipStack[cs] || x > this.clipStack[cs] + this.clipStack[cs + 2])) continue;
      const y = placed ? l.ly[j] : l.y[i];
      this.map(x, y);
      const m = o.nextLetter();
      o.codes[m] = l.codes[j];
      o.x[m] = this.px;
      o.ly[m] = this.py;
      o.angle[m] = (this.tangent ? this.angleAt(x) : 0) + (placed ? l.angle[j] : 0);
      o.vs[m] = placed ? l.vs[j] : 1;
      o.colors[m] = l.color[i] !== NoColor ? NoColor : l.colors[j];
    }
    o.end();
  }

  /** A text at its anchor's place, upright. */
  private text(out: PassWriter, p: PassList, i: number) {
    const l = p.texts;
    this.map(l.x[i], l.y[i]);
    out.pass.op(Op.text);
    const o = out.pass.texts;
    const k = o.next();
    o.str[k] = l.str[i];
    o.x[k] = this.px;
    o.y[k] = this.py;
    o.anchor[k] = l.anchor[i];
    o.font[k] = l.font[i];
    o.color[k] = l.color[i];
    o.width[k] = l.width[i];
    o.maxWidth[k] = l.maxWidth[i];
  }
}
