/**
 * Copyright (c) 2026 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

// Draw lists (_spec/README.md "Render", Draw lists): CSS pixels, origin top left; a list per pass, drawn in order;
// primitives as structs of arrays in growable typed arrays kept frame to frame (a reset only zeroes the counts, so
// nothing is allocated once warm). Within a pass the primitives draw in the order they were added: `ops` holds each
// one's kind, and each kind's arrays are read with a cursor of their own.

/**
 * The passes, below to above: the Canvas2D backend draws the first two on one canvas, the emphasis on an overlay. The
 * `window` pass is the first of the emphasis (a ring's or a navigator's window, under the highlights).
 */
export const Passes = ['background', 'data', 'window', 'highlight', 'selection', 'hover', 'cursor', 'labels'] as const;

export type PassName = (typeof Passes)[number];

/** The passes a hover, a selection, the cursor or a window's move redraws alone. */
export const EmphasisPasses: readonly PassName[] = ['window', 'highlight', 'selection', 'hover', 'cursor'];

export const Op = { rect: 0, path: 1, glyphs: 2, text: 3, image: 4, clip: 5, unclip: 6 } as const;

export type OpKind = keyof typeof Op;

/** Rounded ends of a squircle rect: an end cut by the view or a block is square (the frame rule). */
export const Round = { none: 0, left: 1, right: 2, both: 3 } as const;

export const PathCommand = { move: 0, line: 1, cubic: 2, close: 3 } as const;

/** Points each path command reads. */
export const PathCommandPoints = [1, 1, 3, 0] as const;

export const Anchor = { start: 0, middle: 1, end: 2 } as const;

/** No color: an absent stroke, or a glyph run colored per letter. */
export const NoColor = 0;

/** `a`, or a larger copy of it with its contents kept (buffers grow while they fill). */
export function grown<T extends Float64Array | Float32Array | Int32Array | Uint32Array | Uint16Array | Uint8Array>(
  a: T,
  need: number,
): T {
  if (need <= a.length) return a;
  const next = new (a.constructor as new (n: number) => T)(Math.max(need, a.length * 2));
  next.set(a);
  return next;
}

/** `rect`: x, y, w, h, fill, stroke and its width, a squircle radius (0: square) and which ends are round. */
export class RectList {
  n = 0;
  x = new Float32Array(64);
  y = new Float32Array(64);
  w = new Float32Array(64);
  h = new Float32Array(64);
  fill = new Uint32Array(64);
  stroke = new Uint32Array(64);
  strokeWidth = new Float32Array(64);
  radius = new Float32Array(64);
  round = new Uint8Array(64);

  push(
    x: number,
    y: number,
    w: number,
    h: number,
    fill: number,
    stroke = NoColor,
    strokeWidth = 0,
    radius = 0,
    round = 3,
  ) {
    const i = this.n++;
    if (i >= this.x.length) this.grow(i + 1);
    this.x[i] = x;
    this.y[i] = y;
    this.w[i] = w;
    this.h[i] = h;
    this.fill[i] = fill;
    this.stroke[i] = stroke;
    this.strokeWidth[i] = strokeWidth;
    this.radius[i] = radius;
    this.round[i] = radius > 0 ? round : 0;
    return i;
  }

  /** The next rect's index, its arrays grown to hold it (the writer fills them in). */
  next(): number {
    const i = this.n++;
    if (i >= this.x.length) this.grow(i + 1);
    return i;
  }

  private grow(need: number) {
    this.x = grown(this.x, need);
    this.y = grown(this.y, need);
    this.w = grown(this.w, need);
    this.h = grown(this.h, need);
    this.fill = grown(this.fill, need);
    this.stroke = grown(this.stroke, need);
    this.strokeWidth = grown(this.strokeWidth, need);
    this.radius = grown(this.radius, need);
    this.round = grown(this.round, need);
  }
}

/** `path`: commands and points in shared arrays (move, line, cubic, close), fill, stroke and its width, dash. */
export class PathList {
  n = 0;
  cmdStart = new Uint32Array(16);
  cmdCount = new Uint32Array(16);
  pointStart = new Uint32Array(16);
  fill = new Uint32Array(16);
  stroke = new Uint32Array(16);
  strokeWidth = new Float32Array(16);
  dash = new Float32Array(16);
  commands = new Uint8Array(64);
  commandCount = 0;
  /** x, y pairs. */
  points = new Float32Array(128);
  pointCount = 0;

  begin() {
    const i = this.n;
    if (i >= this.cmdStart.length) {
      const need = i + 1;
      this.cmdStart = grown(this.cmdStart, need);
      this.cmdCount = grown(this.cmdCount, need);
      this.pointStart = grown(this.pointStart, need);
      this.fill = grown(this.fill, need);
      this.stroke = grown(this.stroke, need);
      this.strokeWidth = grown(this.strokeWidth, need);
      this.dash = grown(this.dash, need);
    }
    this.cmdStart[i] = this.commandCount;
    this.pointStart[i] = this.pointCount;
  }

  command(c: number) {
    this.commands = grown(this.commands, this.commandCount + 1);
    this.commands[this.commandCount++] = c;
  }

  point(x: number, y: number) {
    this.points = grown(this.points, 2 * this.pointCount + 2);
    this.points[2 * this.pointCount] = x;
    this.points[2 * this.pointCount + 1] = y;
    this.pointCount++;
  }

  end(fill: number, stroke = NoColor, strokeWidth = 0, dash = 0) {
    const i = this.n++;
    this.cmdCount[i] = this.commandCount - this.cmdStart[i];
    this.fill[i] = fill;
    this.stroke[i] = stroke;
    this.strokeWidth[i] = strokeWidth;
    this.dash[i] = dash;
    return i;
  }
}

/**
 * `glyphs`: a run of single letters in one font, an x (center) per letter, one baseline y, a color per letter or one.
 * A **placed** run (`placed`) gives each letter its own baseline `ly`, angle (radians, clockwise, about its place) and
 * vertical scale `vs` (letters along a ring, a logo's stacked letters, a row's letters).
 */
export class GlyphList {
  n = 0;
  font = new Uint32Array(16);
  y = new Float32Array(16);
  /** The run's color; `NoColor`: each letter's own. */
  color = new Uint32Array(16);
  start = new Uint32Array(16);
  count = new Uint32Array(16);
  /** 1: each letter has its own y, angle and vertical scale. */
  placed = new Uint8Array(16);
  /** Per letter, shared by the runs: its UTF-16 code, center x, color. */
  codes = new Uint16Array(256);
  x = new Float32Array(256);
  colors = new Uint32Array(256);
  /** Per letter of a placed run: its baseline, angle (radians, clockwise) and vertical scale. */
  ly = new Float32Array(256);
  angle = new Float32Array(256);
  vs = new Float32Array(256);
  letters = 0;

  begin(font: number, y: number, color: number) {
    const i = this.n;
    if (i >= this.font.length) this.growRuns(i + 1);
    this.font[i] = font;
    this.y[i] = y;
    this.color[i] = color;
    this.start[i] = this.letters;
    this.count[i] = 0;
    this.placed[i] = 0;
  }

  letter(code: number, x: number, color = NoColor) {
    const j = this.letters++;
    if (j >= this.codes.length) this.growLetters(j + 1);
    this.codes[j] = code;
    this.x[j] = x;
    this.colors[j] = color;
    this.count[this.n]++;
    if (this.placed[this.n] === 1) this.plain(j);
  }

  /** A plain letter of a placed run: on the run's baseline, unturned, unscaled. */
  plain(j: number) {
    this.ly[j] = this.y[this.n];
    this.angle[j] = 0;
    this.vs[j] = 1;
  }

  /** Turns the open run placed: the letters it holds so far become plain placed letters. */
  place() {
    const i = this.n;
    if (this.placed[i] === 1) return;
    this.placed[i] = 1;
    for (let j = this.start[i]; j < this.letters; j++) this.plain(j);
  }

  end() {
    return this.n++;
  }

  /** The next run's index (begun by the writer), its arrays grown to hold it. */
  nextRun(): number {
    const i = this.n;
    if (i >= this.font.length) this.growRuns(i + 1);
    this.start[i] = this.letters;
    this.count[i] = 0;
    this.placed[i] = 0;
    return i;
  }

  /** The next letter's index in the open run, its arrays grown to hold it. */
  nextLetter(): number {
    const j = this.letters++;
    if (j >= this.codes.length) this.growLetters(j + 1);
    this.count[this.n]++;
    // a plain letter in a placed run (a placed letter's numbers are set after)
    if (this.placed[this.n] === 1) this.plain(j);
    return j;
  }

  private growRuns(need: number) {
    this.font = grown(this.font, need);
    this.y = grown(this.y, need);
    this.color = grown(this.color, need);
    this.start = grown(this.start, need);
    this.count = grown(this.count, need);
    this.placed = grown(this.placed, need);
  }

  private growLetters(need: number) {
    this.codes = grown(this.codes, need);
    this.x = grown(this.x, need);
    this.colors = grown(this.colors, need);
    this.ly = grown(this.ly, need);
    this.angle = grown(this.angle, need);
    this.vs = grown(this.vs, need);
  }
}

/** `text`: an interned string at x, y (baseline) with an anchor, a font, a color, its width, a max width (faded past). */
export class TextList {
  n = 0;
  str = new Uint32Array(16);
  x = new Float32Array(16);
  y = new Float32Array(16);
  anchor = new Uint8Array(16);
  font = new Uint32Array(16);
  color = new Uint32Array(16);
  /** Its measured width: a backend fades it out past `maxWidth` without measuring. */
  width = new Float32Array(16);
  /** Infinity: never faded. */
  maxWidth = new Float32Array(16);

  push(
    str: number,
    x: number,
    y: number,
    anchor: number,
    font: number,
    color: number,
    width: number,
    maxWidth: number,
  ) {
    const i = this.n++;
    if (i >= this.str.length) this.grow(i + 1);
    this.str[i] = str;
    this.x[i] = x;
    this.y[i] = y;
    this.anchor[i] = anchor;
    this.font[i] = font;
    this.color[i] = color;
    this.width[i] = width;
    this.maxWidth[i] = maxWidth;
    return i;
  }

  /** The next text's index, its arrays grown to hold it. */
  next(): number {
    const i = this.n++;
    if (i >= this.str.length) this.grow(i + 1);
    return i;
  }

  private grow(need: number) {
    this.str = grown(this.str, need);
    this.x = grown(this.x, need);
    this.y = grown(this.y, need);
    this.anchor = grown(this.anchor, need);
    this.font = grown(this.font, need);
    this.color = grown(this.color, need);
    this.width = grown(this.width, need);
    this.maxWidth = grown(this.maxWidth, need);
  }
}

/** RGBA bytes, row by row: dense strips and matrices. The layer that made it owns and reuses the bytes. */
export interface ImageData {
  readonly width: number;
  readonly height: number;
  readonly data: Uint8Array;
}

/** `image`: placed at x, y, w, h (CSS px), its pixels scaled to the box. */
export class ImageList {
  n = 0;
  x = new Float32Array(8);
  y = new Float32Array(8);
  w = new Float32Array(8);
  h = new Float32Array(8);
  /** The images, `n` of them (a slot past `n` is emptied at a reset: the list holds no image it no longer draws). */
  images: ImageData[] = [];

  /** The next image's index, its arrays grown to hold it. */
  next(): number {
    const i = this.n++;
    if (i >= this.x.length) {
      this.x = grown(this.x, i + 1);
      this.y = grown(this.y, i + 1);
      this.w = grown(this.w, i + 1);
      this.h = grown(this.h, i + 1);
    }
    return i;
  }

  push(x: number, y: number, w: number, h: number, image: ImageData) {
    const i = this.n++;
    if (i >= this.x.length) {
      this.x = grown(this.x, i + 1);
      this.y = grown(this.y, i + 1);
      this.w = grown(this.w, i + 1);
      this.h = grown(this.h, i + 1);
    }
    this.x[i] = x;
    this.y[i] = y;
    this.w[i] = w;
    this.h[i] = h;
    this.images[i] = image;
    return i;
  }
}

/**
 * `clip` pushes a rect (intersected with the clip around it), a squircle where it has a radius (a track's well: its
 * round ends `round` names, the others square); `unclip` pops it.
 */
export class ClipList {
  n = 0;
  x = new Float32Array(16);
  y = new Float32Array(16);
  w = new Float32Array(16);
  h = new Float32Array(16);
  radius = new Float32Array(16);
  round = new Uint8Array(16);

  push(x: number, y: number, w: number, h: number, radius = 0, round: number = Round.both) {
    const i = this.n++;
    if (i >= this.x.length) {
      this.x = grown(this.x, i + 1);
      this.y = grown(this.y, i + 1);
      this.w = grown(this.w, i + 1);
      this.h = grown(this.h, i + 1);
      this.radius = grown(this.radius, i + 1);
      this.round = grown(this.round, i + 1);
    }
    this.x[i] = x;
    this.y[i] = y;
    this.w[i] = w;
    this.h[i] = h;
    this.radius[i] = radius > 0 && round !== 0 ? radius : 0;
    this.round[i] = radius > 0 ? round : 0;
    return i;
  }
}

/** One pass: its primitives in the order added. */
export class PassList {
  ops = new Uint8Array(256);
  opCount = 0;
  readonly rects = new RectList();
  readonly paths = new PathList();
  readonly glyphs = new GlyphList();
  readonly texts = new TextList();
  readonly images = new ImageList();
  readonly clips = new ClipList();
  /** Clips pushed and not popped. */
  depth = 0;

  /** The pass's index in `Passes`. */
  readonly index: number;

  constructor(readonly name: PassName) {
    this.index = Passes.indexOf(name);
  }

  op(kind: number) {
    if (this.opCount >= this.ops.length) this.ops = grown(this.ops, this.opCount + 1);
    this.ops[this.opCount++] = kind;
  }

  /** The primitives in the pass, clip pushes and pops included. */
  get size(): number {
    return this.opCount;
  }

  reset() {
    this.opCount = 0;
    this.depth = 0;
    this.rects.n = 0;
    this.paths.n = 0;
    this.paths.commandCount = 0;
    this.paths.pointCount = 0;
    this.glyphs.n = 0;
    this.glyphs.letters = 0;
    this.texts.n = 0;
    // the images let go, their slots kept (a length of 0 would let the array's storage go, made again next frame)
    const images = this.images.images as (ImageData | undefined)[];
    for (let i = 0; i < this.images.n; i++) images[i] = undefined;
    this.images.n = 0;
    this.clips.n = 0;
  }
}

/** Reads a pass's primitives in order: each callback gets its list and the primitive's index there. */
export interface PassVisitor {
  rect(l: RectList, i: number): void;
  path(l: PathList, i: number): void;
  glyphs(l: GlyphList, i: number): void;
  text(l: TextList, i: number): void;
  image(l: ImageList, i: number): void;
  clip(l: ClipList, i: number): void;
  unclip(): void;
}

export function walkPass(p: PassList, v: PassVisitor): void {
  let r = 0;
  let pa = 0;
  let g = 0;
  let t = 0;
  let im = 0;
  let c = 0;
  for (let k = 0; k < p.opCount; k++) {
    switch (p.ops[k]) {
      case Op.rect:
        v.rect(p.rects, r++);
        break;
      case Op.path:
        v.path(p.paths, pa++);
        break;
      case Op.glyphs:
        v.glyphs(p.glyphs, g++);
        break;
      case Op.text:
        v.text(p.texts, t++);
        break;
      case Op.image:
        v.image(p.images, im++);
        break;
      case Op.clip:
        v.clip(p.clips, c++);
        break;
      case Op.unclip:
        v.unclip();
        break;
    }
  }
}

/** The interned strings draw lists keep across frames by default (`RenderConfig.internedStrings`). */
export const DefaultMaxStrings = 4096;

/**
 * The draw lists of a frame: a list per pass, the interned strings (texts and fonts) and the size they were laid out
 * for. Reused frame to frame: `reset` keeps every buffer. `key` is what they were built from (the renderer's).
 */
export class DrawLists {
  readonly passes: Readonly<Record<PassName, PassList>>;
  width = 0;
  height = 0;
  key = '';
  private readonly _strings: string[] = [];
  private readonly _ids = new Map<string, number>();

  /** @param maxStrings the interned strings kept across frames (`RenderConfig.internedStrings`). */
  constructor(readonly maxStrings = DefaultMaxStrings) {
    const passes = {} as Record<PassName, PassList>;
    for (const name of Passes) passes[name] = new PassList(name);
    this.passes = passes;
  }

  /**
   * Clears every pass (or the ones named) for a new frame. The interned strings stay (a frame drawn again interns
   * nothing new, so it allocates nothing); a full reset clears them once they are more than `maxStrings`.
   */
  reset(passes: readonly PassName[] = Passes): void {
    for (let i = 0; i < passes.length; i++) this.passes[passes[i]].reset();
    if (passes.length === Passes.length && this._strings.length > this.maxStrings) {
      this._strings.length = 0;
      this._ids.clear();
    }
  }

  /** The id of an interned string (a text, a font). */
  intern(s: string): number {
    const id = this._ids.get(s);
    return id === undefined ? this.add(s) : id;
  }

  private add(s: string): number {
    const id = this._strings.length;
    this._strings.push(s);
    this._ids.set(s, id);
    return id;
  }

  string(id: number): string {
    return this._strings[id];
  }

  private readonly _writers = new Map<PassName, PassWriter>();

  /** The writer of a pass: one a pass, kept, so drawing allocates none. */
  writer(pass: PassName): PassWriter {
    let w = this._writers.get(pass);
    if (!w) this._writers.set(pass, (w = new PassWriter(this, this.passes[pass])));
    return w;
  }

  /** Primitives in every pass. */
  get size(): number {
    let n = 0;
    for (const p of Passes) n += this.passes[p].size;
    return n;
  }
}

/** Adds primitives to one pass of draw lists: what layer types draw through. */
export class PassWriter {
  /**
   * A rect's, an image's or a text's numbers, set by the caller before `rectArgs`, `imageArgs` or `textArgs`: x, y, w,
   * h, fill, stroke, stroke width, radius, round ends (a text's: x, y, width, max width, color). A call boxes each number it is passed where the compiler did not inline it (an
   * allocation); an array's numbers are never boxed, so the emphasis passes, drawn rarely next to the data passes and
   * so never sure to be inlined, draw through it (nothing allocated once warm).
   */
  readonly args = new Float64Array(9);

  constructor(
    readonly lists: DrawLists,
    readonly pass: PassList,
  ) {}

  /** The rect `args` holds (all nine numbers). */
  rectArgs() {
    const a = this.args;
    if (a[2] <= 0 || a[3] <= 0) return;
    this.pass.op(Op.rect);
    const l = this.pass.rects;
    const i = l.next();
    l.x[i] = a[0];
    l.y[i] = a[1];
    l.w[i] = a[2];
    l.h[i] = a[3];
    l.fill[i] = a[4];
    l.stroke[i] = a[5];
    l.strokeWidth[i] = a[6];
    l.radius[i] = a[7];
    l.round[i] = a[7] > 0 ? a[8] : 0;
  }

  /** Begins a glyph run at the baseline `args[1]` holds (`beginGlyphs`, its number through `args`). */
  beginGlyphsArgs(font: string, color = NoColor) {
    this.pass.op(Op.glyphs);
    const l = this.pass.glyphs;
    const i = l.nextRun();
    l.font[i] = this.lists.intern(font);
    l.y[i] = this.args[1];
    l.color[i] = color;
  }

  /**
   * A text at the numbers `args` holds: x, y (baseline), its measured width, its max width (faded past it), its color
   * (a color past a small integer passed to a call is boxed: it goes through `args` too).
   */
  textArgs(s: string, anchor: number, font: string) {
    const a = this.args;
    this.pass.op(Op.text);
    const l = this.pass.texts;
    const i = l.next();
    l.str[i] = this.lists.intern(s);
    l.x[i] = a[0];
    l.y[i] = a[1];
    l.anchor[i] = anchor;
    l.font[i] = this.lists.intern(font);
    l.color[i] = a[4];
    l.width[i] = a[2];
    l.maxWidth[i] = a[3];
  }

  /** The image at the box `args` holds (its first four numbers). */
  imageArgs(image: ImageData) {
    const a = this.args;
    this.pass.op(Op.image);
    const l = this.pass.images;
    const i = l.next();
    l.x[i] = a[0];
    l.y[i] = a[1];
    l.w[i] = a[2];
    l.h[i] = a[3];
    l.images[i] = image;
  }

  rect(
    x: number,
    y: number,
    w: number,
    h: number,
    fill: number,
    stroke = NoColor,
    strokeWidth = 0,
    radius = 0,
    round = 3,
  ) {
    if (w <= 0 || h <= 0) return;
    this.pass.op(Op.rect);
    // stored here, not passed on: a number passed to a call that is not inlined is boxed (nothing allocated once warm)
    const l = this.pass.rects;
    const i = l.next();
    l.x[i] = x;
    l.y[i] = y;
    l.w[i] = w;
    l.h[i] = h;
    l.fill[i] = fill;
    l.stroke[i] = stroke;
    l.strokeWidth[i] = strokeWidth;
    l.radius[i] = radius;
    l.round[i] = radius > 0 ? round : 0;
  }

  /**
   * A filled square rect, no stroke: `rect` with its defaults, small, so the compiler inlines it into the loops that
   * call it (a call that is not inlined boxes the numbers it is passed).
   */
  fill(x: number, y: number, w: number, h: number, fill: number) {
    if (w <= 0 || h <= 0) return;
    this.pass.op(Op.rect);
    const l = this.pass.rects;
    const i = l.next();
    l.x[i] = x;
    l.y[i] = y;
    l.w[i] = w;
    l.h[i] = h;
    l.fill[i] = fill;
    l.stroke[i] = NoColor;
    l.strokeWidth[i] = 0;
    l.radius[i] = 0;
    l.round[i] = 0;
  }

  /** A path: `build` adds its commands through `move`, `line`, `cubic`, `close`. */
  path(build: (p: PathList) => void, fill: number, stroke = NoColor, strokeWidth = 0, dash = 0) {
    this.pass.op(Op.path);
    const p = this.pass.paths;
    p.begin();
    build(p);
    p.end(fill, stroke, strokeWidth, dash);
  }

  /**
   * A path through the `n` points of `xy` (x, y pairs): a move to the first, lines to the rest, closed when `close`; its
   * fill, stroke, stroke width and dash are `args[4]` to `args[7]` (a number passed to a call that is not inlined is
   * boxed: the data passes of the small layers draw through here, nothing allocated once warm).
   */
  pathArgs(xy: Float64Array, n: number, close: boolean) {
    if (n < 2) return;
    const a = this.args;
    this.pass.op(Op.path);
    const p = this.pass.paths;
    p.begin();
    p.commands = grown(p.commands, p.commandCount + n + 1);
    p.points = grown(p.points, 2 * (p.pointCount + n));
    for (let i = 0; i < n; i++) {
      p.commands[p.commandCount++] = i === 0 ? PathCommand.move : PathCommand.line;
      p.points[2 * p.pointCount] = xy[2 * i];
      p.points[2 * p.pointCount + 1] = xy[2 * i + 1];
      p.pointCount++;
    }
    if (close) p.commands[p.commandCount++] = PathCommand.close;
    const k = p.n++;
    p.cmdCount[k] = p.commandCount - p.cmdStart[k];
    p.fill[k] = a[4];
    p.stroke[k] = a[5];
    p.strokeWidth[k] = a[6];
    p.dash[k] = a[7];
  }

  /**
   * A cubic curve from `xy[o], xy[o + 1]` through the controls at `o + 2` and `o + 4` to `o + 6`, stroked: its fill,
   * stroke, stroke width and dash are `args[4]` to `args[7]` (an arc of pairs; nothing allocated once warm).
   */
  curveArgs(xy: Float64Array, o: number) {
    const a = this.args;
    this.pass.op(Op.path);
    const p = this.pass.paths;
    p.begin();
    p.commands = grown(p.commands, p.commandCount + 2);
    p.points = grown(p.points, 2 * (p.pointCount + 4));
    p.commands[p.commandCount++] = PathCommand.move;
    p.commands[p.commandCount++] = PathCommand.cubic;
    for (let i = 0; i < 4; i++) {
      p.points[2 * p.pointCount] = xy[o + 2 * i];
      p.points[2 * p.pointCount + 1] = xy[o + 2 * i + 1];
      p.pointCount++;
    }
    const k = p.n++;
    p.cmdCount[k] = p.commandCount - p.cmdStart[k];
    p.fill[k] = a[4];
    p.stroke[k] = a[5];
    p.strokeWidth[k] = a[6];
    p.dash[k] = a[7];
  }

  /** Begins a glyph run; `letter` adds to it, `endGlyphs` closes it. */
  beginGlyphs(font: string, y: number, color = NoColor) {
    this.pass.op(Op.glyphs);
    const l = this.pass.glyphs;
    const i = l.nextRun();
    l.font[i] = this.lists.intern(font);
    l.y[i] = y;
    l.color[i] = color;
  }

  /** A letter at the x `args[0]` holds, in the color `args[4]` holds (`letter`, its numbers through `args`). */
  letterArgs(code: number) {
    const l = this.pass.glyphs;
    const j = l.nextLetter();
    l.codes[j] = code;
    l.x[j] = this.args[0];
    l.colors[j] = this.args[4];
  }

  letter(code: number, x: number, color = NoColor) {
    const l = this.pass.glyphs;
    const j = l.nextLetter();
    l.codes[j] = code;
    l.x[j] = x;
    l.colors[j] = color;
  }

  /**
   * A placed letter (the open run becomes placed): x, y, angle (radians, clockwise), vertical scale and color from
   * `args[0]` to `args[4]` (its numbers through `args`, as `letterArgs`).
   */
  placedLetterArgs(code: number) {
    const l = this.pass.glyphs;
    l.place();
    const j = l.nextLetter();
    const a = this.args;
    l.codes[j] = code;
    l.x[j] = a[0];
    l.ly[j] = a[1];
    l.angle[j] = a[2];
    l.vs[j] = a[3];
    l.colors[j] = a[4];
  }

  endGlyphs() {
    this.pass.glyphs.end();
  }

  text(
    s: string,
    x: number,
    y: number,
    anchor: number,
    font: string,
    color: number,
    width: number,
    maxWidth = Infinity,
  ) {
    this.pass.op(Op.text);
    const l = this.pass.texts;
    const i = l.next();
    l.str[i] = this.lists.intern(s);
    l.x[i] = x;
    l.y[i] = y;
    l.anchor[i] = anchor;
    l.font[i] = this.lists.intern(font);
    l.color[i] = color;
    l.width[i] = width;
    l.maxWidth[i] = maxWidth;
  }

  image(x: number, y: number, w: number, h: number, image: ImageData) {
    this.pass.op(Op.image);
    this.pass.images.push(x, y, w, h, image);
  }

  clip(x: number, y: number, w: number, h: number, radius = 0, round: number = Round.both) {
    this.pass.op(Op.clip);
    this.pass.clips.push(x, y, w, h, radius, round);
    this.pass.depth++;
  }

  /** A clip to the box `args` holds (its first four numbers): fractional numbers passed to a call are boxed. */
  clipArgs() {
    const a = this.args;
    this.pass.op(Op.clip);
    const l = this.pass.clips;
    if (l.n >= l.x.length) l.push(a[0], a[1], a[2], a[3]);
    else {
      const i = l.n++;
      l.x[i] = a[0];
      l.y[i] = a[1];
      l.w[i] = a[2];
      l.h[i] = a[3];
      l.radius[i] = 0;
      l.round[i] = 0;
    }
    this.pass.depth++;
  }

  /**
   * A clip to the squircle `args` holds (x, y, w, h, its radius, its round ends: a track's well), a rect where the
   * radius is 0 or no end is round.
   */
  clipSquircleArgs() {
    const a = this.args;
    this.pass.op(Op.clip);
    const l = this.pass.clips;
    if (l.n >= l.x.length) l.push(a[0], a[1], a[2], a[3], a[4], a[5]);
    else {
      const i = l.n++;
      l.x[i] = a[0];
      l.y[i] = a[1];
      l.w[i] = a[2];
      l.h[i] = a[3];
      const r = a[4] > 0 && a[5] !== 0 ? a[4] : 0;
      l.radius[i] = r;
      l.round[i] = r > 0 ? a[5] : 0;
    }
    this.pass.depth++;
  }

  unclip() {
    if (this.pass.depth === 0) return;
    this.pass.op(Op.unclip);
    this.pass.depth--;
  }
}

/** A pass's draw function of a layout (a layer type's, or the view's own). */
export type PassDraw<L, E> = (lay: L, out: PassWriter, emphasis: E) => void;

/** A layout's draw functions by pass, at each pass's index in `Passes` (`passDraws`); none for a pass it does not draw. */
export type PassDraws<L, E> = readonly (PassDraw<L, E> | undefined)[];

/** Draw functions by pass name, as a table by pass index (`drawPass` reads it). */
export function passDraws<L, E>(byPass: Partial<Record<PassName, PassDraw<L, E>>>): PassDraws<L, E> {
  return Passes.map((p) => byPass[p]);
}

/**
 * Draws the writer's pass of a layout with its function in the table (_spec/README.md "Render", Draw lists: nothing
 * allocated once warm). The function is read from an array and called at this one site, which every pass of the view
 * and of every layer type goes through: no call target is known to the compiler here, so no pass function is inlined
 * into its caller, and each is compiled on its own with the whole inlining budget for the writer's calls inside it. A
 * writer call left out of its caller's inlining boxes each number it is passed: an allocation a primitive.
 */
export function drawPass<L, E>(draws: PassDraws<L, E>, lay: L, out: PassWriter, emphasis: E): void {
  const draw = draws[out.pass.index];
  if (draw !== undefined) draw(lay, out, emphasis);
}

/** Path command helpers (on a `PathList` being built). */
export function moveTo(p: PathList, x: number, y: number) {
  p.command(PathCommand.move);
  p.point(x, y);
}

export function lineTo(p: PathList, x: number, y: number) {
  p.command(PathCommand.line);
  p.point(x, y);
}

export function cubicTo(p: PathList, x1: number, y1: number, x2: number, y2: number, x: number, y: number) {
  p.command(PathCommand.cubic);
  p.point(x1, y1);
  p.point(x2, y2);
  p.point(x, y);
}

export function closePath(p: PathList) {
  p.command(PathCommand.close);
}
