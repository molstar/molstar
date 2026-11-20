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
// their angles, a one-row image (a density strip) its runs of one color as sectors, a taller image refused (its layout
// says so), a clip applied in the strip's frame before warping. Buffers are kept: nothing allocated once warm.

import { grown, NoColor, Op, PassList, Passes, PassWriter, PathCommand, type PathList } from '../draw-lists.js';
import { fontSize } from '../config.js';
import { squirclePathArgs } from '../squircle.js';

/** A letter's center above its baseline, in ems (the layers' `fontSize * 0.35` from a cell's middle to its baseline). */
const CenterEm = 0.35;

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

  /**
   * Moves what was drawn since `open` by strip `i`'s translation (`dx[i]`, `dy[i]`): read here, never passed (a
   * fractional number passed to a call that is not inlined is boxed: a text's line tops are fractional).
   */
  close(out: PassWriter, dxs: Float64Array, dys: Float64Array, i: number): void {
    const dx = dxs[i];
    const dy = dys[i];
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
    for (let j = this.clips; j < c.n; j++) {
      c.x[j] += dx;
      c.y[j] += dy;
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
  /** Primitives refused (an image of more than one row): its layout reports them. */
  refused = 0;
  /** A scratch pass a pass, named as it (what a strip draws reads its pass's name: the emphasis does). */
  private readonly scratches = Passes.map((name) => new PassList(name));
  private readonly writers: (PassWriter | undefined)[] = [];
  private scratch = this.scratches[0];
  /** The clip rects open in the strip's frame (x, y, w, h each). */
  private clipStack = new Float64Array(16);
  private clipDepth = 0;
  /** The square's outline a rect is traced through, and the rect's numbers. */
  private readonly square = new PassList('data');
  /** A one-row image's runs, a rect each, warped as rects are. */
  private readonly runs = new PassList('data');
  private readonly rectArgs = new Float64Array(6);

  /** The page point of a frame point: into `px`, `py`. */
  px = 0;
  py = 0;
  map(x: number, y: number): void {
    this.qx = x;
    this.qy = y;
    this.mapQ();
  }

  /** Appends (`px`, `py`) to a path's points: no number passed (a call that is not inlined boxes its numbers). */
  private put(o: PathList): void {
    const k = o.pointCount;
    o.points = grown(o.points, 2 * k + 2);
    o.points[2 * k] = this.px;
    o.points[2 * k + 1] = this.py;
    o.pointCount = k + 1;
  }

  /** The frame point the draw loops map, set as fields: a number passed to a call that is not inlined is boxed. */
  private qx = 0;
  private qy = 0;
  /** `(qx, qy)` to the page, into `px`, `py`, and its angle into `qt`. */
  private qt = 0;
  private mapQ(): void {
    const t = this.a0 + ((this.qx - this.x0) / this.width) * this.sweep;
    const r = this.r0 - this.qy;
    this.qt = t;
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
          if (!this.band(out, p, im++)) this.refused++;
          break;
        case Op.clip:
          // the clip by its index: a number passed to a call that is not inlined is boxed
          this.pushClip(p, c++);
          break;
        case Op.unclip:
          if (this.clipDepth > 0) this.clipDepth--;
          break;
      }
    }
    p.reset();
  }

  private pushClip(p: PassList, c: number) {
    const l = p.clips;
    let x0 = l.x[c];
    let y0 = l.y[c];
    let x1 = x0 + l.w[c];
    let y1 = y0 + l.h[c];
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

  /**
   * A one-row image (a density strip: a color a device column) as the ring band it stands for, each run of one color an
   * annular sector, so a ring draws a layer's density as a track does; false for a taller image (a matrix).
   */
  private band(out: PassWriter, p: PassList, i: number): boolean {
    const l = p.images;
    const image = l.images[i];
    if (!image || image.height !== 1) return false;
    const d = image.data;
    const n = image.width;
    const cw = l.w[i] / Math.max(1, n);
    const runs = this.runs;
    for (let c = 0; c < n; ) {
      const k = 4 * c;
      if (d[k + 3] === 0) {
        c++;
        continue;
      }
      let e = c + 1;
      while (
        e < n &&
        d[4 * e] === d[k] &&
        d[4 * e + 1] === d[k + 1] &&
        d[4 * e + 2] === d[k + 2] &&
        d[4 * e + 3] === d[k + 3]
      )
        e++;
      runs.reset();
      runs.rects.push(
        l.x[i] + c * cw,
        l.y[i],
        (e - c) * cw,
        l.h[i],
        ((d[k] << 24) | (d[k + 1] << 16) | (d[k + 2] << 8) | d[k + 3]) >>> 0,
        NoColor,
        0,
        0,
        0,
      );
      this.rect(out, runs, 0);
      c = e;
    }
    return true;
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
    const qp = q.paths;
    qp.begin();
    if (l.round[i] === 0 || l.radius[i] <= 0) {
      // a square rect (most marks at a ring's density) traced here, its corners through `px`, `py`: nothing boxed
      qp.command(PathCommand.move);
      this.px = x;
      this.py = y;
      this.put(qp);
      qp.command(PathCommand.line);
      this.px = x + w;
      this.put(qp);
      qp.command(PathCommand.line);
      this.py = y + h;
      this.put(qp);
      qp.command(PathCommand.line);
      this.px = x;
      this.put(qp);
      // back to the start as a line (so it bends), then closed: the squircle's outline with no corners
      qp.command(PathCommand.line);
      this.py = y;
      this.put(qp);
      qp.command(PathCommand.close);
    } else {
      // through an array: a number passed to a call that is not inlined is boxed (nothing allocated once warm)
      const a = this.rectArgs;
      a[0] = x;
      a[1] = y;
      a[2] = w;
      a[3] = h;
      a[4] = l.radius[i];
      a[5] = l.round[i];
      squirclePathArgs(qp, a);
    }
    // the style written in place: a color passed to a call is boxed
    const qk = q.paths.end(NoColor);
    q.paths.fill[qk] = l.fill[i];
    q.paths.stroke[qk] = l.stroke[i];
    q.paths.strokeWidth[qk] = l.strokeWidth[i];
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
        this.qx = lx;
        this.qy = ly;
        this.mapQ();
        o.command(PathCommand.move);
        this.put(o);
      } else if (cmd === PathCommand.line) {
        const x = l.points[2 * pt];
        const y = l.points[2 * pt + 1];
        pt++;
        // the segment's length along the ring at its radius, cut every `step` px so it bends
        const r = Math.max(0, this.r0 - (ly + y) / 2);
        const arc = Math.abs(((x - lx) / this.width) * this.sweep) * r;
        const n = Math.max(1, Math.min(512, Math.ceil(arc / this.step)));
        for (let k = 1; k <= n; k++) {
          this.qx = lx + ((x - lx) * k) / n;
          this.qy = ly + ((y - ly) * k) / n;
          this.mapQ();
          o.command(PathCommand.line);
          this.put(o);
        }
        lx = x;
        ly = y;
      } else if (cmd === PathCommand.cubic) {
        o.command(PathCommand.cubic);
        for (let k = 0; k < 3; k++, pt++) {
          this.qx = l.points[2 * pt];
          this.qy = l.points[2 * pt + 1];
          this.mapQ();
          this.put(o);
        }
        lx = l.points[2 * pt - 2];
        ly = l.points[2 * pt - 1];
      } else o.command(PathCommand.close);
    }
    const k = o.end(NoColor);
    o.fill[k] = l.fill[i];
    o.stroke[k] = l.stroke[i];
    o.strokeWidth[k] = l.strokeWidth[i];
    o.dash[k] = l.dash[i];
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
    const size = this.sizeOf(out, l.font[i]);
    // a letter is kept while any of it is in the clip: half a mono advance (0.3 em) either side of its center
    const half = cs >= 0 ? 0.3 * size : 0;
    // its center is above its baseline by `CenterEm` (the layers set a letter's baseline at its cell's middle plus it)
    const up = CenterEm * size;
    for (let j = l.start[i], end = j + l.count[i]; j < end; j++) {
      const x = l.x[j];
      if (cs >= 0 && (x + half <= this.clipStack[cs] || x - half >= this.clipStack[cs] + this.clipStack[cs + 2]))
        continue;
      const a = placed ? l.angle[j] : 0;
      const v = placed ? l.vs[j] : 1;
      // the letter's center is what is warped (its cell's angle, its band's radius), then its baseline is put back
      // below it in its own turned frame: a baseline warped instead puts the letter's center off its cell, along the
      // ring where the letter stands upright
      this.qx = x + Math.sin(a) * v * up;
      this.qy = (placed ? l.ly[j] : l.y[i]) - Math.cos(a) * v * up;
      this.mapQ();
      const t = (this.tangent ? this.qt : 0) + a;
      const m = o.nextLetter();
      o.codes[m] = l.codes[j];
      o.x[m] = this.px - Math.sin(t) * v * up;
      o.ly[m] = this.py + Math.cos(t) * v * up;
      o.angle[m] = t;
      o.vs[m] = v;
      o.colors[m] = l.color[i] !== NoColor ? NoColor : l.colors[j];
    }
    o.end();
  }

  /** A font's size in px by its interned id, read once a font. */
  private readonly sizes = new Map<number, number>();
  private sizeOf(out: PassWriter, font: number): number {
    let v = this.sizes.get(font);
    if (v === undefined) this.sizes.set(font, (v = fontSize(out.lists.string(font))));
    return v;
  }

  /** A text at its anchor's place, upright. */
  private text(out: PassWriter, p: PassList, i: number) {
    const l = p.texts;
    this.qx = l.x[i];
    this.qy = l.y[i];
    this.mapQ();
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
